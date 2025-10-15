
# -*- coding: utf-8 -*-
"""
AISLAMIENTO CRIOGÉNICO - PUR (Opción 2: Script Python ejecutable)
Autor: ChatGPT (para Fernando Peluffo)
Licencia: MIT

Descripción
-----------
Calcula, por diámetro de columna/tubería, el perfil de temperatura radial, el flujo de calor
y el peso del sistema de aislamiento PUR. Usa solo pulgadas para LONGITUD radial.
Salidas clave por cada caso:
- q (W/m) y (BTU/h·ft)
- q'' en la superficie externa (W/m²)
- Perfil T(r) a través del sistema (gráfico con fondo blanco, líneas gruesas y marcadores numerados)
- Peso del PUR por metro (kg/m) y por longitud especificada (kg)

Cómo usar
---------
1) Edite "config_aislamiento.json" en la misma carpeta para ingresar sus casos.
2) Ejecute:
   > python aislamiento_pur.py
   (Necesita Python 3 + numpy + matplotlib; si no tiene Python, puede indicarme y le preparo un .exe)
3) Los resultados se guardan en ./salidas/:
   - PNG por cada diámetro (gráfica de temperatura radial)
   - CSV resumen de resultados
   - TXT de consola (bitácora)

Unidades y convenciones
-----------------------
- Distancia radial y diámetros: pulgadas [in]
- Longitud axial para cálculo de q: metros [m] (interno), pero se reporta q también en BTU/h·ft
- Temperaturas: °C (entrada/salida). Si quiere °F, ajuste "mostrar_temp_F" en el config.
- Convección interna/externa: h en W/m²·K
- Conductividad PUR dependiente de la temperatura (k(T)) aproximada linealmente: k = a + b*(T-25)
  con a=0.022 W/m·K y b=1.5e-4 W/m·K/°C (ajustable en el config).

Limitaciones
------------
- Modelo 1D, régimen permanente, k(PUR) efectiva evaluada iterativamente a T_media en PUR.
- Pared metálica interna opcional (suma de resistencia por conducción cilíndrica con k_acero).
- Barreras/vainas adicionales pueden incluirse como capas extra (ver "capas_extra").
"""

import json
import math
import os
from dataclasses import dataclass, asdict
from typing import List, Tuple

import numpy as np
import matplotlib.pyplot as plt


# ----------------------------- Utilidades de conversión -----------------------------
IN_TO_M = 0.0254
FT_TO_M = 0.3048
W_TO_BTU_PER_H = 3.412141633  # 1 W = 3.412 BTU/h

def btu_per_h_per_ft_from_W_per_m(q_W_per_m: float) -> float:
    # q [W/m] -> [BTU/h·ft]
    return q_W_per_m * W_TO_BTU_PER_H / (1/FT_TO_M)

# ----------------------------- Datos y modelos -----------------------------
@dataclass
class ConfigCaso:
    nombre: str
    diametro_interno_in: float          # diámetro interno "fluido" [in]
    espesor_pared_in: float = 0.0       # espesor pared metálica [in] (opcional)
    k_acero_W_mK: float = 16.0          # conductividad del acero [W/m·K]
    T_fluido_C: float = -160.0          # °C
    T_ambiente_C: float = 25.0          # °C
    h_in_W_m2K: float = 500.0           # convección interna [W/m²·K] (ebullición/forzada)
    h_out_W_m2K: float = 10.0           # convección externa [W/m²·K] (aire quieto ≈ 5-15)
    espesor_PUR_in: float = 2.0         # espesor PUR [in]
    densidad_PUR_kg_m3: float = 35.0    # kg/m³
    largo_calculo_m: float = 1.0        # longitud axial para q [m]
    mostrar_temp_F: bool = False        # además de °C
    # k(T) PUR = a + b*(T-25°C)
    kPUR_a_W_mK: float = 0.022
    kPUR_b_W_mK_perC: float = 1.5e-4
    # Capas extra: lista de dicts con {"nombre":str,"espesor_in":float,"k_W_mK":float}
    capas_extra: list = None

@dataclass
class ResultadoCaso:
    nombre: str
    q_W_m: float
    q_BTU_h_ft: float
    qpp_W_m2: float
    peso_PUR_kg_m: float
    peso_PUR_kg_total: float
    r_in_in: float
    r_ext_PUR_in: float

def k_PUR_func(T_C: float, a: float, b: float) -> float:
    return a + b * (T_C - 25.0)

def resistencia_conveccion(h_W_m2K: float, r_m: float, L_m: float) -> float:
    # R_conv (cilíndrica) ~ 1/(h*A) con A=2π r L
    A = 2.0 * math.pi * r_m * L_m
    return 1.0 / (h_W_m2K * A) if h_W_m2K > 0 else 0.0

def resistencia_conduccion_cil(k_W_mK: float, r1_m: float, r2_m: float, L_m: float) -> float:
    return math.log(r2_m / r1_m) / (2.0 * math.pi * k_W_mK * L_m)

def perfil_T_cilindro(q_W_m: float, nodos_r_m: np.ndarray, r_interfaces_m: List[Tuple[float, float, float]],
                      T_in_C: float, R_conv_in: float) -> np.ndarray:
    """
    Construye perfil T(r) desde r_in hasta r_out, restando caída de T por cada resistencia.
    r_interfaces_m: lista de (r1, r2, k) para regiones sólidas en orden radial.
    """
    T_r = np.zeros_like(nodos_r_m)
    T_actual = T_in_C - q_W_m * R_conv_in  # caída por convección interna
    r_prev = nodos_r_m[0]
    T_r[0] = T_actual
    # Para cada segmento sólido, computar caída local proporcional a ln(r/r1)
    idx = 1
    r1_global = nodos_r_m[0]
    for (r1, r2, k) in r_interfaces_m:
        # puntos dentro de este segmento
        while idx < len(nodos_r_m) and nodos_r_m[idx] <= r2 + 1e-12:
            r = nodos_r_m[idx]
            frac = math.log(r / r1) / math.log(r2 / r1) if r > r1 else 0.0
            deltaT_seg = q_W_m * resistencia_conduccion_cil(k, r1, r2, 1.0)  # L=1 para deltaT relativo
            T_r[idx] = T_actual - frac * deltaT_seg
            idx += 1
        # Al final del segmento, actualizamos T_actual
        T_actual = T_r[idx-1] if idx-1 >=0 else T_actual
    # Los nodos restantes (si quedaran por redondeo) se rellenan con T_actual
    while idx < len(nodos_r_m):
        T_r[idx] = T_actual
        idx += 1
    return T_r

def resolver_caso(cfg: ConfigCaso) -> Tuple[ResultadoCaso, np.ndarray, np.ndarray]:
    # Radios en m
    r_in_m = (cfg.diametro_interno_in * IN_TO_M) / 2.0
    r_pared_ext_m = r_in_m + (cfg.espesor_pared_in * IN_TO_M)
    r_PUR_ext_m = r_pared_ext_m + (cfg.espesor_PUR_in * IN_TO_M)

    # Resistencias
    R_conv_in = resistencia_conveccion(cfg.h_in_W_m2K, r_in_m, cfg.largo_calculo_m)

    R_pared = 0.0
    if cfg.espesor_pared_in > 0 and cfg.k_acero_W_mK > 0:
        R_pared = resistencia_conduccion_cil(cfg.k_acero_W_mK, r_in_m, r_pared_ext_m, cfg.largo_calculo_m)

    # Capas extra (conducción)
    capas = []
    if cfg.capas_extra:
        r1 = r_PUR_ext_m  # suponemos se agregan después del PUR
        for capa in cfg.capas_extra:
            r2 = r1 + capa["espesor_in"] * IN_TO_M
            capas.append((r1, r2, capa["k_W_mK"]))
            r1 = r2

    # Iteración para k_PUR efectiva (dependiente de T)
    T_in = cfg.T_fluido_C
    T_out = cfg.T_ambiente_C

    # Inicial: k_PUR a T media estimada
    T_media_PUR = (T_in + T_out) / 2.0
    kPUR = k_PUR_func(T_media_PUR, cfg.kPUR_a_W_mK, cfg.kPUR_b_W_mK_perC)

    # Construir R_total en función de kPUR
    def R_total(kpur: float) -> float:
        R_PUR = resistencia_conduccion_cil(kpur, r_pared_ext_m, r_PUR_ext_m, cfg.largo_calculo_m)
        R_capas = sum(resistencia_conduccion_cil(k, r1, r2, cfg.largo_calculo_m) for (r1, r2, k) in capas)
        R_conv_out = resistencia_conveccion(cfg.h_out_W_m2K, (capas[-1][1] if capas else r_PUR_ext_m), cfg.largo_calculo_m)
        return R_conv_in + R_pared + R_PUR + R_capas + R_conv_out

    # Iteración simple de punto fijo para k(T_media_PUR)
    for _ in range(20):
        Rtot = R_total(kPUR)
        q = (T_in - T_out) / Rtot
        # temperatura en la pared externa del PUR (antes de capas extra)
        R_hasta_PURext = R_conv_in + R_pared + resistencia_conduccion_cil(kPUR, r_pared_ext_m, r_PUR_ext_m, cfg.largo_calculo_m)
        T_PUR_ext = T_in - q * R_hasta_PURext
        T_media_PUR = (T_in + T_PUR_ext) / 2.0
        kPUR_new = k_PUR_func(T_media_PUR, cfg.kPUR_a_W_mK, cfg.kPUR_b_W_mK_perC)
        if abs(kPUR_new - kPUR) < 1e-5:
            kPUR = kPUR_new
            break
        kPUR = kPUR_new

    # Resultados finales
    Rtot = R_total(kPUR)
    q_W_m = (T_in - T_out) / Rtot
    q_BTU_h_ft = btu_per_h_per_ft_from_W_per_m(q_W_m)

    # Flujo en la superficie exterior
    r_ext_real_m = capas[-1][1] if capas else r_PUR_ext_m
    A_ext = 2.0 * math.pi * r_ext_real_m * cfg.largo_calculo_m
    # Temperatura en la superficie externa sólida antes de convección externa
    R_hasta_extsolida = Rtot - resistencia_conveccion(cfg.h_out_W_m2K, r_ext_real_m, cfg.largo_calculo_m)
    T_ext_solida = T_in - q_W_m * R_hasta_extsolida
    qpp_W_m2 = cfg.h_out_W_m2K * (T_ext_solida - T_out)

    # Peso del PUR
    area_anular_PUR_m2 = math.pi * (r_PUR_ext_m**2 - r_pared_ext_m**2)
    peso_PUR_kg_m = cfg.densidad_PUR_kg_m3 * area_anular_PUR_m2
    peso_PUR_kg_total = peso_PUR_kg_m * cfg.largo_calculo_m

    # Construir perfil radial
    r0 = r_in_m
    r_list = [r0]
    etiquetas_segmentos = []
    interfaces = []  # (r1, r2, k) para perfil_T
    # pared metálica (opcional)
    if cfg.espesor_pared_in > 0:
        interfaces.append((r_in_m, r_pared_ext_m, cfg.k_acero_W_mK))
        r_list += list(np.linspace(r_in_m, r_pared_ext_m, 20)[1:])
        etiquetas_segmentos.append("Pared metálica")
        r0 = r_pared_ext_m
    # PUR
    interfaces.append((r0, r_PUR_ext_m, kPUR))
    r_list += list(np.linspace(r0, r_PUR_ext_m, 50)[1:])
    etiquetas_segmentos.append("PUR (k≈{:.4f} W/m·K)".format(kPUR))
    r0 = r_PUR_ext_m
    # Capas extra
    for (r1, r2, k) in capas:
        interfaces.append((r1, r2, k))
        r_list += list(np.linspace(r1, r2, 20)[1:])
        etiquetas_segmentos.append("Capa extra k={:.3f}".format(k))
        r0 = r2
    nodos_r_m = np.array(r_list)
    T_r = np.zeros_like(nodos_r_m)
    # Perfil T: caídas en cada tramo + convección interna inicial
    T_r = perfil_T_cilindro(q_W_m, nodos_r_m, interfaces, T_in, R_conv_in)

    # Preparar ejes en pulgadas y temperatura °C (y opcional °F)
    r_in_in = cfg.diametro_interno_in / 2.0
    r_ext_PUR_in = r_in_in + cfg.espesor_pared_in + cfg.espesor_PUR_in
    r_in_pulg = nodos_r_m / IN_TO_M
    T_C = T_r
    # Gráfica
    fig = plt.figure(figsize=(7.5, 5.0))
    plt.plot(r_in_pulg, T_C, linewidth=3.0, marker='o', markersize=3)
    # Marcadores numerados cada ~10 puntos para legibilidad
    for i, (x, y) in enumerate(zip(r_in_pulg[::10], T_C[::10]), start=1):
        plt.annotate(str(i), (x, y), textcoords="offset points", xytext=(5,5), fontsize=8)
    plt.title(f"SECCIÓN TÍPICA AISLAMIENTO CRIOGÉNICO - PUR\n{cfg.nombre}", fontsize=11)
    plt.xlabel("Radio [in]")
    ylabel = "Temperatura [°C]"
    plt.ylabel(ylabel)
    plt.grid(True, linewidth=0.4, alpha=0.4)
    plt.tight_layout()
    plt.gcf().patch.set_facecolor('white')
    out_dir = os.path.join("salidas")
    os.makedirs(out_dir, exist_ok=True)
    png_path = os.path.join(out_dir, f"{cfg.nombre}_perfil_T.png")
    plt.savefig(png_path, dpi=160, bbox_inches="tight")
    plt.close(fig)

    res = ResultadoCaso(
        nombre=cfg.nombre,
        q_W_m=q_W_m,
        q_BTU_h_ft=q_BTU_h_ft,
        qpp_W_m2=qpp_W_m2,
        peso_PUR_kg_m=peso_PUR_kg_m,
        peso_PUR_kg_total=peso_PUR_kg_total,
        r_in_in=r_in_in,
        r_ext_PUR_in=r_ext_PUR_in
    )
    # Guardar TXT resumen del caso
    txt_path = os.path.join(out_dir, f"{cfg.nombre}_resumen.txt")
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write(f"=== {cfg.nombre} ===\n")
        f.write(f"q = {q_W_m:.2f} W/m  |  {q_BTU_h_ft:.2f} BTU/h·ft\n")
        f.write(f"q''_ext = {qpp_W_m2:.2f} W/m²\n")
        f.write(f"Peso PUR = {peso_PUR_kg_m:.3f} kg/m  |  en {cfg.largo_calculo_m:.2f} m: {peso_PUR_kg_total:.3f} kg\n")
        f.write(f"Radio interno = {r_in_in:.3f} in  |  Radio externo PUR = {r_ext_PUR_in:.3f} in\n")
        f.write(f"Imagen perfil: {png_path}\n")

    return res, r_in_pulg, T_C


def cargar_config(path_json: str) -> List[ConfigCaso]:
    with open(path_json, "r", encoding="utf-8") as f:
        data = json.load(f)
    casos = []
    for item in data["casos"]:
        cfg = ConfigCaso(**item)
        casos.append(cfg)
    return casos

def main():
    print("== AISLAMIENTO CRIOGÉNICO - PUR | Opción 2 (Python) ==")
    print("Leyendo configuración: config_aislamiento.json")
    casos = cargar_config("config_aislamiento.json")
    out_dir = os.path.join("salidas")
    os.makedirs(out_dir, exist_ok=True)
    resumen_csv = os.path.join(out_dir, "resumen_resultados.csv")
    with open(resumen_csv, "w", encoding="utf-8") as fcsv:
        fcsv.write("nombre,q_W_m,q_BTU_h_ft,qpp_W_m2,peso_PUR_kg_m,peso_PUR_kg_total,r_in_in,r_ext_PUR_in\n")
        for cfg in casos:
            res, r_in_pulg, T_C = resolver_caso(cfg)
            fcsv.write("{},{:.6f},{:.6f},{:.6f},{:.6f},{:.6f},{:.6f},{:.6f}\n".format(
                res.nombre, res.q_W_m, res.q_BTU_h_ft, res.qpp_W_m2,
                res.peso_PUR_kg_m, res.peso_PUR_kg_total, res.r_in_in, res.r_ext_PUR_in
            ))
            print(f"[OK] {res.nombre}: q={res.q_W_m:.2f} W/m ({res.q_BTU_h_ft:.2f} BTU/h·ft), Peso PUR={res.peso_PUR_kg_m:.3f} kg/m")
    print(f"\nListo. Revise la carpeta '{out_dir}'.")

if __name__ == "__main__":
    main()

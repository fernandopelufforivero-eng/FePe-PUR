AISLAMIENTO CRIOGÉNICO – SISTEMA PUR | FePe – Cryogenic PUR System
Documento técnico de despliegue y operación web
Autor: F. Peluffo – Octubre 2025
Revisión: 01

1.0 OBJETO DEL DOCUMENTO
Este documento define el propósito técnico del sistema web “FePe – Cryogenic PUR System”, herramienta para diseño de espesor de aislamiento PUR en equipos cilíndricos criogénicos. El módulo calcula espesor total, número de capas y ubicación de barrera de vapor, y genera la lámina técnica para impresión en Carta (8.5×11 in).

Nota técnica (principio físico):
Conducción radial estacionaria (Ley de Fourier) en cilindros, con resistencias en serie (convección interna, acero, PUR y convección externa). El flujo lineal es:
Q = 2πL (Ti − To) / [ 1/(hi ri) + ln(r2/r1)/k_acero + ln(r3/r2)/k_PUR + 1/(ho r3) ]
El modelo ajusta k_PUR = a + b (T − 25°C) iterando a temperatura media. Resultados en W/m, BTU/h·ft y W/m².

2.0 CONTENIDO DEL PAQUETE
- index.html           Aplicación web autosuficiente (HTML5/JS) con modo DISEÑO (espesor, capas y BV).
- README_DEPLOY.txt    Guía técnica de despliegue (GitHub Pages, Netlify, dominio propio).

3.0 PROCEDIMIENTO DE DESPLIEGUE
3.1 GitHub Pages
  a) Crear repositorio (p.ej. FePe-PUR) y subir index.html.
  b) Settings → Pages → Deploy from branch → main / root. Guardar.
  c) Abrir la URL publicada en Safari/Chrome.
3.2 Netlify
  a) https://app.netlify.com/drop → soltar index.html.
  b) Netlify publica un enlace del tipo: https://NOMBRE.netlify.app
3.3 Dominio propio
  a) Subir index.html al directorio público (/public_html o raíz del hosting).
  b) Acceder por https://su-dominio/...

4.0 CONDICIONES DE OPERACIÓN Y COMPATIBILIDAD
- Compatible con Safari (iOS/iPadOS), Chrome, Edge y Firefox.
- Formato de impresión: Carta (8.5×11 in), fondo blanco, firma automática “F. Peluffo – Octubre 2025”.
- Requisitos de navegador: JavaScript (ES6) y Canvas 2D.
- No requiere Python ni instalaciones locales.

5.0 CONTROL DE VERSIONES Y FIRMA TÉCNICA
F. Peluffo – Octubre 2025
Revisión: 01

© 2025 FePe – Cryogenic PUR System. Todos los derechos reservados.

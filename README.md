# Scripts-pipeline-SVs
Protocolo bioinformático para el procesado de los archivos VCF de 3 "variant callers" de variantes estructurales en muestras de secuenciación de genoma completo (WGS) y su posterior puesta en común de las variantes detectadas y recogida en un único archivo BED. 

# ¿Como correr esta pipeline?
Los archivos de la muestra que se desee analizar debe de estar en un mismo directorio organizados de la siguiente manera:
> [!Por ejemplo]
- /0032_0032/
    - Manta/
          Manta_0032.vcf
    - delly/
          delly_0032.bcf
    - GRIDSS/
          GRIDSS_0032.vcf

De esta forma solo se tiene que ejecutar el siguiente comando:
./pipe_SV.sh ruta/directorio/0032_0032
> [!IMPORTANT]
> la ruta al directorio debe de escribirse sin el "/" final, de lo contrario dará problemas.

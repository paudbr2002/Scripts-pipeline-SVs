# Scripts-pipeline-SVs
Protocolo bioinformático para el procesado de los archivos VCF de 3 "variant callers" de variantes estructurales en muestras de secuenciación de genoma completo (WGS) y su posterior puesta en común de las variantes detectadas y recogida en un único archivo BED. 

## ¿Como correr esta pipeline?
Los archivos de la muestra que se desee analizar debe de estar en un mismo directorio organizados de la siguiente manera:
[Por ejemplo]
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
> la ruta al directorio se debe escribir sin el "/" final.

## Esquema funcionamiento general ##

![. Flujo de trabajo de los scripts incorporados para procesar los VCF de cada programa y obtener un único archivo con
todas las variantes incluyendo las detectadas por varios programas. Los scripts destinados al procesado de los VCF se recuadran en
rojo y los destinados a la fusión de las variantes en un solo archivo en naranja. DEL = Deleción, INS= Inserción, DUP =
Duplicación, TRA = Translocación, INV= Inversión.
](https://github.com/paudbr2002/Scripts-pipeline-SVs/blob/main/Esquema_pipeline.png)
Inicialmente se crean directorios para filtrar las variantes con calidad “PASS”, se transforman los archivos VCF en BED y vía diferentes scripts se obtienen los archivos formateados de cada programa para cada SV.
Posteriormente se separan los SV por genotipo y se ordenan y agrupan las variantes por coordenadas. Se utilizan scripts específicos para procesar diferentes tipos de SV (en inversiones y deleciones se
utiliza la librería bedtools merge) y se añaden las anotaciones necesarias. Finalmente, todos los archivos se combinan en un archivo final con un formato específico (Tabla Suplementaria 1), que contiene todas las variantes agrupadas por coordenadas y la información sobre el solapamiento de la variante de cada programa con el consenso.
## Versiones ##
bcftools 1.15, bedtools 2.27.1, BEDOPS v2.4.41, Python 3.12.3, Bash 5.0.17, R 4.4.0
## Autora ##
Paula de Blas Rioja 2024

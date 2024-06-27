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
](https://www.canva.com/design/DAGJUFlLyz8/9X5Crutbigm7ziRjYaAePQ/edit)

## Versiones ##
bcftools 1.15, bedtools 2.27.1, BEDOPS v2.4.41, Python 3.12.3, Bash 5.0.17, R 4.4.0

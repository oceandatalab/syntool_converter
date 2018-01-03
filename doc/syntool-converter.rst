syntool-converter
=================
The aim of the syntool-converter command is to handle the specificities of the
input data (units conversion, sensor correction, etc…) and convert them in
GeoTIFF or NetCDF so that the other Syntool softwares only have to support
these two formats.

::
    syntool-converter -t TYPE -i INPUT -o OUTPUT [-h] [-l] [-d DATE]
                      [--keep-output][-opt OPTIONS] [--options-file OPT_FILE]|


+---------------------------------+----------------------+-----------------------------------------------------------+
| Parameter                       | Format               | Description                                               |
+=================================+=======================+==========================================================+
| -i INPUT, --input INPUT         | path                  | Read data from the INPUT file or directory (depends on   |
|                                 |                       | the data reader).                                        |
+---------------------------------+-----------------------+----------------------------------------------------------+
| -o OUTPUT, --output OUTPUT      | path                  | Save results in the OUTPUT directory (will be created if |
|                                 |                       | needed).                                                 |
+---------------------------------+-----------------------+----------------------------------------------------------+
| -t TYPE                         | string                | Use a reader which supports TYPE data as input. Complete |
|                                 |                       | list of supported types can be obtained by displaying    |
|                                 |                       | help with the -h option.                                 |
+---------------------------------+-----------------------+----------------------------------------------------------+
| -opt OPTIONS, --options OPTIONS | key=value [key=value] | TYPE-specific option that will be passed to the          |
|                                 |                       | conversion method. OPTIONS is a list of key-value        |
|                                 |                       | couples (format is “key=value”) separated by a blank     |
|                                 |                       | space.                                                   |
+---------------------------------+-----------------------+----------------------------------------------------------+
| --options-file OPT_FILE         | path                  | Replaces the -t and -opt values passed using the command |
|                                 |                       | line by the ones contained in the OPT_FILE text file.    |
+---------------------------------+-----------------------+----------------------------------------------------------+
| -d DATE, --date DATE            | YYYYmmddTHHMMSS       | Applicable only for input files that contain several     |
|                                 |                       | granules: extract data for the granule which corresponds |
|                                 |                       | to the specified DATE.                                   |
+---------------------------------+-----------------------+----------------------------------------------------------+
| -l, --list                      |                       | Parse INPUT as a text file where each line is to be      |
|                                 |                       | treated as the path of a input file.                     |
+---------------------------------+-----------------------+----------------------------------------------------------+
| --keep-output                   |                       | Save results directly in OUTPUT (no subdirectory named   |
|                                 |                       | after the product).                                      |
+---------------------------------+-----------------------+----------------------------------------------------------+                                                      
| -h, --help                      |                       | Display help message.                                    |
+---------------------------------+-----------------------+----------------------------------------------------------+  

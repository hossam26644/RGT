# RGT
Repeats Genotyping Tool


```
RGT
│─── main.py               
│─── RGT.py                
│─── settings.json    
│─── README.md
│ 
└──────────── interface
│             │───  interface.py
│             │───  check_files.py
│             │───  json_parser.py
│             └───  __init__.py
│ 
└──────────── FileReader
│             │─── ReadFile.py 
│             └───  __init__.py
│
└──────────── ExcelExporter
│             │───  ExcelExport.py 
│             └───  __init__.py
│ 
└──────────── genotyper
│             │───  GenoType.py
│             │───  Repeat.py
│             │───  SmartString.py
│             │───  GroupingString.py
│             │───  revComplementry.py
│             └───  __init__.py
│ 
└──────────── graphsplotter
│             │───  plotter.py
│             │───  table_2d_plotter.py   
│             │───  plot_3D.py
│             └───  __init__.py
│ 
└──────────── allelesdetector
              │───  AllelesDetector.py
              │───  MatchingSequence.py
              │───  PeakIdentifier.py
              └───  __init__.py

```

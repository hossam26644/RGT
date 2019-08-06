# RGT
Repeats Genotyping Tool

## This is a developers guide, if you are looking for a user manual: sorry i did not make one yet

RGT is a software to:
1. Extract strucutre from fastq files representing targeted sequencing samples of SSRs (simple sequence repeats)
2. Export the identified germline alleles of the sample (assuming it is a human sample)

## The approach taken here is an alignment free approach
* Because targeted sequencing gives thousands of reads, why do alignment
* We just extract structures and get the most repeating one
* This of course means you are loosing a percentage of the reads as any sequencing error will result in a mistakenly identified structure, but who cares if you have tens of thousands of reads

## Strucutres are identified by the expected repeating units identified by the user
* strucutres are stored with their number of repeating units
* A graph of the number of repeating units count and the abundance of samples having the same count is plotted
* Peaks from this graph is used to identified the germline alleles

> Every structure identified is stored regardless of the sequencing errors and variations and exported in a separate excel file for the user

##### The most abundant repeat strucutres are matched with the peaks identified in the "repeating units count vs structures having this count plot" to get the two human germline alleles

## lets get to the developers guide:
#### first the strucutre of the software:
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
└──────────── filereader
│             │─── ReadFile.py 
│             └───  __init__.py
│
└──────────── excelexporter
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
## how parallel processing is done:
* The main class uses joblib parallel, it does the job for you
* Processes are created by it and instances of the RGT class where the analysis is done are excuted in parallel
* The results of all processes are put in a list of dictionaries, there is a method *get_collective_dictionary_from_list_of_output_dictionaries* (do you like my naming style?) to (as you would imagine) get a collective dictionary from the list of output dictionaries

## Do I have to explain every single detail? sorry I am not going to
* the interface gets user inputs, nothing smart in that 
* excelexporter exports tables to excel sheets, nothing smart also I guess
* file reader parses fastq files and extracts the reads:
  * it searches for the flanking sequences and extract reads from between these flanks
  * flanking sequence match is identified when a sequence has a [hamming distance](https://en.wikipedia.org/wiki/Hamming_distance) of one or zero with the flank
  * when there is no upstream flank found (if uyser specified a start flank) the whole read is discarded
  * if the user specidfies to discard reads when there is no down stream flank.. and no down stream flank is found: guess what, the read is discarded
  * The percentage of the discarded reads is calculated and the sample is flagged if this percentage is greater than a threshold identified by the user (or default value set by the software)

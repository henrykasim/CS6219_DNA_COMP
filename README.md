# CS6219_DNA_COMP
## Objective:
Primer Library Design. on given primers length, 'CG' content limit, and edit distance among the primers limit


## Prerequisite:

``` pip install python-levenshtein ```


## Usage:
``` python main.py ```

## Usage (with parameters):
``` python main.py -m 1 -l 20 -min 0.45 -max 0.55 -d 0.4 ```

* -m 	is method to generate primers library (default is 1)
* -l    is strands length (default is 8)
* -min  is percentage of strands length as the minimum number of character 'C' and 'G' (default is 0.45)
* -max  is percentage of strands length as the maximum number of character 'C' and 'G' (default is 0.55)
* -d    is percentage of strands length as the minimum edit distance allowed (default is 0.4)

## Available Method:
* Method 1: Generate all possible strands and filter based on CG_limit and ED_limit (edit distance). 

* Method 2: Generate selected CG_allowed strands and evaluate strands based on least similar sum of Edit distance

* Method 3: Genetic Algorithm with Cross over, mutation, shift, shuffle, and remove some population.



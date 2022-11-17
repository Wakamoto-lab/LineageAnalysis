# What is this
This repository provides MATLAB codes used in "A unified framework for measuring selection on cellular lineages and traits" by Shunpei Yamauchi, Takashi Nozoe, Reiko Okura, Edo Kussell and Yuichi Wakaoto

eLife2022;11:e72299 DOI: https://doi.org/10.7554/eLife.72299

# Requirement
- Matlab 2021a
- Matlab Statistics and Machine Learning Toolbox

# Usage
The three files in the cumulant folder correspond to the analysis for each figure as follows;

- exponential/matlab: Figure 5 and its suppemental figures.
  - run "main.m" and chose a directory containing data to be analysed.   
- rpoS/matlab: Figure 6.
  - run "data_preparation.m" to align the format of the Lineage Data, including the number of slices, and then run "analysis.m" to evaluate the quantities in the lineage statistics. 
- stationary/matlab: Figure 7 and its supplemental figures.
  - run "main.m" and chose a directory containing data to be analysed.   

# Author
- Shunpei Yamauchi
- shunpei_yamauchi@cell.c.u-tokyo.ac.jp

# License
This is under [MIT license](https://en.wikipedia.org/wiki/MIT_License).

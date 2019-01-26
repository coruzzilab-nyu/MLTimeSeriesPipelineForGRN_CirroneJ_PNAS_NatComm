Tutorial Instructions ValidateFinalRanking_and_PRCurve_Step

- install anaconda: https://conda.io/miniconda.html

- then from the command line, you have to create a conda environment as follows:
python2 -m pip install ipykernel
python2 -m ipykernel install --user
conda create -n ipykernel_py2 python=2 ipykernel
source activate ipykernel_py2    # On Windows, remove the word 'source'
python -m ipykernel install --user
pip install numpy, pandas, matplotlib, sklearn

- Finally type "jupyter notebook": a browser page will open up.

- Open "ValidateFinalRanking_And_PRCurve.ipynb"

- "#Set the most important parameters: 
pvalue_thres(float), list_name(string), remove_noExpInProtopl(boolean), remove_noJITs(boolean) and then Run the jupyter cell with shift_enter"
( you just have to set those params in the cell accordingly and run the notebook cell with "shift+enter")


- Copy to "TFTargets folder" the list of files related to the TFs you want to validate.
(You can find an example inside "TFTargets folder", please copy files with the same format as the example)

- "results" folder will contain the output results

- PlotPackage_For_Final_AUPR_Figure for plotting the final PR_curve for the paper.

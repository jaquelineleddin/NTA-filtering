NTAWorkflow is a R package for filtering feature list exported from MS-DIAL after processing LC-HRMS data in non-target analysis.

---

## Instructions

### Installing R and RStudio
If you don't have R and RStudio on your PC, then you need to install it. 
1. Go to https://cran.uni-muenster.de and follow the instructions to install R
2. Go to https://posit.co/download/rstudio-desktop/ to install RStudio

### Install Packages
1. Start RStudio.
2. Create a new project. Go to the menu bar: File -> New Project -> Existing Directory. Browse for the NTAWorkflow directory. When you go on *Files* on the right side in RStudio you will see all files of the directory. We will come to those files later.
3. Go to *Packages* on the right side. Click on *install*, then choose *Install from: Package Archive File*. Then browse for the file path, where you have saved the NTAWorkflow.tar.gz file. After installing, you will find NTAWorkflow in the list of packages.
4. You will need some packages from the internet. In RStudio go to *Packages* on the right side. Click on *Install* and choose *Install from: Repository (CRAN)*. For *Packages* type in: ggplot2, gridExtra, shiny. Then click on *install*. After installing you will see ggplot2, gridExtra, shiny in the list in RStudio under *Packages*

---

Now there are two ways of using the scripts: with the notebook or with the app. For beginners it is easier to use the app.

### a) Filtering data with the app

1. Start RStudio, if it's not started yet.
2. On the right side click on *Files*. If you don't see the files from the NTAWorkflow directory, go to *File* in the menu bar on top of the window and *Open project* and browse for the NTAWorkflow directory on your PC and click on the file .RProj.
3. In *Files* on the right side click on *NTA App* and open ui.R
4. The script will be open in the top left window. Don't write anything in this file! In this window on the top right side appeared the button *Run App*. Click on it to start the app.

### b) Filtering data with the R Notebook

For every dataset you should create a directory and put your exported feature list of MS-DIAL in there. Copy the downlowded notebook (notebook_NTA_filtered.Rmd) into this directory. 
Start RStudio, go to: File -> new Project on the top left. Choose *Existing Directory*, browse the directory with the notebook, then go on *Create Project*. If you click on the right side *Files* you should see the notebook and your raw data text file. Double click on the notebook will open it. Then follow the instructions for the filtering process in the directory **NTANotebook**.



# 2011-jet-inclusivecrosssection-analysis

Only works with ROOT6

Before you start

  1. Create empty directory for output plots
  2. Checkout and compile the tool RooUnfold which will help you unfold the data 

```
$ mkdir plots
$ cd code
$ svn co https://svnsrv.desy.de/public/unfolding/RooUnfold/trunk RooUnfold
$ cd RooUnfold
$ make
$ cd ..

```

Run the full analysis. If you want to run the code with data, change nothing. If you want to run the code with MC, open the file [settings.h](code/settings.h) and change 

```
std::string _jp_type = "DATA"
```
to
```
std::string _jp_type = "MC"
``` 
After, while in the directory [code](code) run

```
$ root -l -b -q mk_runFullAnalysis.C
```

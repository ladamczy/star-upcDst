# star-upcDst

![Bilby Stampede](https://cds.cern.ch/record/2288105/files/fig1.png)

## Overview:

*star-upcDst*, is a new framework mainly developed by Jarda Adam, to simplify analysis that are related to forward and UPC physics. One can use this framework in 2 major steps, or only the 2nd step if the upcDst files are generated by the 1st step already:

- Step 1: Use MuDst files of one's interest as input, run so-called upcDstMaker to produce a upcDst file with information that are needed for one's analysis only. 
- Step 2: Use the result of step 1, the upcDst, to build your personal analysis tree or histograms. 

- *Step 3* (optional): Setup condor job for running large samples

## Step 1:

- Make a clean area on RACF and checkout the repository:

<pre><code> git clone https://github.com/adamjaro/star-upcDst.git </pre></code>

- Go to the main directory star-upcDst:

<pre><code> cd star-upcDst </pre></code>

- Setup the StRoot and build ( already include compiling ) by doing:

<pre><code> ./build-maker.sh </pre></code>

- Next few steps are setting a directory for your input test filelist (below the example is called "test.list"), modify the run config for your need, and then run it:

<pre><code> mkdir txt </pre></code>
<pre><code> cp [your_filelist] ./txt/test.list </pre></code>

- Modify script "config_MC_local.in", there are a few changes needed for your own interest. 1, Modify the "top" directory where your output would go. 2, Modify "add_input", where the second entry is the name you want to give under your "top" directory, and then the path to your test file list. 3. Modify trigger IDs around line 24, where the 2nd and 3rd arguments can be left blank for all run numbers. 4, Modify "is_mc" to be either 0,1,2 for "data","starsim MC", or "embedding MC", respectively.

<pre><code> top path-to-your-own-output-directory </pre></code>

<pre><code> add_input output_folder_name  path-to-filelist </pre></code>

<pre><code> nfiles  3 </pre></code>

<pre><code> add_trigger triggerid </pre></code>

<pre><code> outfile  test.root </pre></code>

and 

<pre><code> is_mc 0 </pre></code>

- Run it to produce the upcDst file

<pre><code> SubmitPlugin.py config_MC_local.in </pre></code>

This will only run interactively for testing purposes

## Step 2:

- With "test.root", the upcDst file, as the input, one can perform simple analysis or run over all events to save information into a tree or histograms. To do that, one needs to compile the package and setup the link to the library. This can be simply done in a few steps:

Under the main directory (under star-upcDst), create a new directory,

<pre><code> mkdir build </pre></code>
<pre><code> cd build </code></pre>
<pre><code> cmake ../ </code></pre>
<pre><code> make </code></pre>

Done. To try an example, go to folder "examples":

<pre><code> cd ../examples </code></pre>

and remember to change the input file name to what has been just created, "test.root", for macro "make_pT.C". Then do, 

<pre><code> root -l run_make_pT.C </code></pre>


- For another example, it uses a C++ code (with ROOT libraries) to read the picoDst file and save output into a small ROOT tree or anything one wants to save. Instead of running it as a .C ROOT macro, this example build an executable and can run as a standard C++ code without using ROOT.

Again, under the main directory, 

<pre><code> cd ./examples/dstreader </code></pre>

make another directory for build,

<pre><code> mkdir build </code></pre>

before one compiles the code, the input file needs to be modified accordingly in "./examples/dstreader/src/AnalysisTemplate.cxx". Again, change it to the "test.root" that has just been created, then go back to directory for build,


<pre><code> cd build </code></pre>
<pre><code> cmake ../ </code></pre>
<pre><code> make </code></pre>

The executable will show up as "AnalysisTemplate", now run it:

<pre><code> ./AnalysisTemplate </code></pre>

It produces an output, "output.root". 

Now it is time for analysis!

## Step 3 (optional): 

There are two different ways to setup your condor jobs. Let's see two examples of submitting jobs for running *Step 1*. Remember, one has to make sure the run script has been setup correctly in *Step 1* if one uses this to produce upcDst files. 

### Example 1:
- MuDst files are on distributed disks and STAR catalog can find your dataset, then one can use Jarda's example as a template. Jarda's example is based on AuAu Run14 UPC events:

This is now done similarly as running interactively, but modifying another config file, e.g., config_JpsiB_run14.in. 
In this config file, please pay attention that it doesn't have the "submit local", meaning that it will submit the jobs on condor. Modify parameters similarly as step 1. 

Finally submit it, 

<pre><code> SubmitPlugin.py config_JpsiB_run14.in </pre></code>

One can check the job status (you can just do condor_q or) using "PrintStat.py":

<pre><code> PrintStat.py -c config_JpsiB_run14.in </pre></code>

Or one can check the job status and also resubmit failed jobs:

<pre><code> PrintStat.py -r -c config_JpsiB_run14.in </pre></code>

PS: -r has to be first then -c.

Finally, to clean up after all outputs are produced correctly, type:

<pre><code> ./clean.sh </pre></code>


### Example 2:
- MuDst files are transfered by someone privately, so it can use a filelist as input to submit condor jobs, then one can use Kong's example as a template. Kong's example is based on dAu Run16 UPC events:

Modify the same config with different "add_input", can be either a filelist or catalog. 


# Complete production guide to star-upcDst

Instructions to run UPC+RP picoDst (upcDst) production are given here,
along with a list of all production options.

## Start the production

Production options are set in a file like  production.in

Command to start the production is:

<pre><code> ./SubmitPlugin.py production.in </pre></code>

executed from  star-upcDst/scheduler/

All the outputs and scheduler files are written to a directory
structure in 'top' directory as set in the options. Outputs
for each input are in the 'top' under the name of a given
input as set with 'add_input' option.

Comments in the file with production options start with #

## Merging the outputs

Output files to each individual input are merged to files
per about 1 GB in size to speed up the analysis and to make
it easier to handle the output.

Using the same file with production options, like the above
example of  production.in,  the merging is started as

<pre><code> ./RunMerge.py production.in </pre></code>

again executed from  star-upcDst/scheduler/

Merged outputs are written to directory specified by 'outdir'
production option, have name derived from 'outfile' option
and list of all produced outputs is written into 'outlist'.

## Overview of running production and resubmission of failed jobs

An overview of ongoing production can be printed as:

<pre><code> ./PrintStat.py -c production.in </pre></code>

To resubmit jobs in error state run the above with '-r'
option:

<pre><code> ./PrintStat.py -r -c production.in </pre></code>

## Example configurations on github

- config_JpsiB_run14.in

Data on J/psi -> e+e- with UPCJpsiB trigger in Run14 AuAu

- config_main_JpsiB_run14.in

The same data on J/psi -> e+e- but also with UPC-main trigger

- config_MC_local.in

Local running with MC files from coherent J/psi embedding

- config_RP_run17.in

Production with Roman Pot data

## Production options

### Required

- top `<top name>`

    main directory for production outputs

- add_input `<input name>` `<catalog query or filelist>`

    inputs to the production, could be catalog query
    or filelist (full path for filelist and should end with .list)

    Could be any number of  add_input  in a single production

- add_trigger <trigger ID> <first run> <last run>

    trigger ID and run range for each ID

    Could be any number of individual trigger IDs. Both first run
    and last run belong to the range.

    The IDs are not considered with MC (both starsim and embedding)

- outfile <name>

    name for merged outputs

- outdir <name>

    directory for merged outputs, relative to top

- outlist <name>

    list of output files after merging


## Contact:
Adam, Jaroslav: <jadam@bnl.gov>









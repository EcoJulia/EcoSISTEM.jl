# SCRC Software checklist 

### Do we have sufficient confidence in the **correctness** of the software to trust the results? 

- Yes / Yes with caveats / No 

> Yes with caveats
>
> The code needs testing against ODEs with both numeric and pen-and-paper solutions in its deterministic mode to give me more confidence.
>
> In stochastic mode, sometimes the tests fail but it's hard to tell if its because of some stochastic variation.

This is your overall judgement on the level of confidence based on all the aspects of the checklist. There is no formulaic way to arrive at this overall assessment based on the individual checklist answers but please explain how the measures in place combine to reach this level of confidence and make clear any caveats (eg applies for certain ways of using the software and not others). 

## Checklist 

For each question: 

*Please explain the situation and include any relevant links (eg tool dashboards, documentation). The sub bullet points are to make the scope of the question clear and should be covered if relevant but do not have to be answered individually.* 

*Please give a “traffic light” assessment: Green: sufficiently addressed, Orange: some outstanding work or caveats about use, Red: Needs to be addressed before software should be relied upon* 

- Can a run be repeated and **reproduce** exactly the same results? 
  - How is stochasticity (randomness) handled? 
  - Is sufficient meta-data logged to enable a run to be reproduced: Exact code version (+ whether the repository was "clean"), versions of libraries (eg.an environment.yml file or similar), all command line arguments, content of any configuration files? 
  - Is there up-to-date documentation which explains precisely how to run the code to reproduce existing results? 

> The code can be repeated exactly by running with the same seed and the same number of threads. It can also be run for deterministic models reproducibly on any number of threads. The seed can be set & recorded by hand. Source code and dependency versions can be recorded in a Project.toml and Manifest.toml file on an ad-hoc basis. No command line arguments are passed into the run scripts for configuration.

- Are there appropriate **tests**? (and are they automated?) 
  - Are there unit tests? What is covered? 
  - System and integration tests? Automated model validation tests? 
  - Regression tests? (showing if changes to the code give unexpected changes in output) 
  - Is there continuous integration? 
  - Is everything you need to run the tests (including documentation) in the repository (or the data pipeline where appropriate)? 

>  The epidemiology parts of the code are well tested with ~70% test coverage. The example run scripts act as end-to-end tests. There is a plan to validate this two dimensional against an ODE code / pen & solutions. There is no regression test as these are easily created ad-hoc during development, and then destroyed. There is continuous integration for linux, macos and windows, which run the tests and deploy the docs. Everything needed to run the tests in in the repo.

- Are the scientific results of runs **robust to different ways of running** the code? *We don't require bitwise identical results here, but the differences should be understandable and not affect the scientific conclusions.*
  - Running on a different machine? 
  - With different number of processes? 
  - With different compilers and optimisation levels? 
  - Running in debug mode? 

> Not tested with lower or higher optimisation levels and/or with debug flags. In stochastic mode, runs with the same seed and number of threads give the same answer respecting the fact that each thread has it’s own random number generator by design. In deterministic mode, the same answer is obtained regardless of the number of threads.

- Has any sort of **automated code checking** been applied? *For C++, this might just be the compiler output when run with "all warnings". It could also be more extensive static analysis. For other languages, it could be e.g. pylint, StaticLint.jl, etc.* 
  - If there are possible issues reported by such a tool, have they all been either fixed or understood not be important? 

> It has not been linted, and no other static analysis tools have been used.

- Is the code clean, generally **understandable and readable** and written according to good software engineering principles? 
  - Is it modular? Are the internal implementation details of one module hidden from other modules? 
  - Commented where necessary? 
  - Avoids red flags such as very long functions, global variables, copy pasted code? 

> The source code is split into three subdirectories: Biodiversity; Climatepref; and Epidemiology. Epideomiology can be considered a dependency of the other two directories (as of today). Some contents of Epidemiology mirror that of the Biodiversity directory but the contents are not unified. Encapsulation could be better. There is room for more code-reuse. Heavy use of arrays and indexing could be improved by with a named index package and/or encapsulation. There are no global variables and, largely, inputs are typed where ambiguous. The intention of outputs are clear. 

- Is there sufficient **documentation**? 
  - Is there a **readme**? 
  - Does the code have **user documentation**? 
  - Does the code have **developer documentation**? 
  - Does the code have **algorithm documentation**? e.g. something that describes how the model is actually simulated, or inference is performed? 
  - Is all the documentation up to date? 

> There is a Readme and user & algorithms docs. No developer docs. There is no documentation per se on how to run the code, but the examples follow Julian style and are self-contained run scripts. 

- Is there suitable **collaboration infrastructure**? 
  - Is the code in a version-controlled repository? 
  - Is there a license? 
  - Is an issue tracker used? 
  - Are there contribution guidelines and review processes? 

> GitHub for version control, issue tracking, travis & appveyor for CI, Codecov, licence included, Copyright (c) 2016-2020: Claire Harris, Richard Reeve and the Scottish COVID-19 Response Consortium. Contributing.md exists.

- Are software **dependencies** listed and of appropriate quality? 
  - Is input and output **data** handled properly? 
  - Does the code **use the data pipeline** for all inputs and outputs? 
  - Is the code **appropriately parameterized** (i.e. have hard coded parameters been removed)? 

> The source code is split into three sub directories: Biodiversity; ClimatePref; and Epidemiology. Epidemiology can be considered a dependency of the other two directories (as of today). Some contents of Epidemiology mirror that of the Biodiversity directory but the contents are not unified. Encapsulation could be better. There is room for more code-reuse. Heavy use of arrays and indexing could be improved by with a named index package and/or encapsulation. There are no global variables and, largely, inputs are typed where ambiguous. The intention of outputs are clear. 

*Detailed* *aspects of data used are covered by the Data checklist and* **not** *here.* 

- Has the model been validated? Is it robust to small changes in input parameters?  

> Not yet, see above.

*These questions are covered by the Model Validation checklist and* **not** *here.* 

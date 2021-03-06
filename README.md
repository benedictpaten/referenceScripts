#Reference Paper Codebase

This is a set of evaluation scripts developed by [Ngan Nguyen ](https://github.com/ngannguyen/) and [Benedict Paten](https://github.com/benedictpaten/) to assess the Reference MHC genome project data.

##Dependencies

* [jobTree](https://github.com/benedictpaten/jobTree)
* [sonLib](https://github.com/benedictpaten/sonLib)
* [cactus](https://github.com/benedictpaten/cactus)
* [assemblaLib](https://github.com/benedictpaten/assemblaLib)

##Running Tests
Download the data set [here](http://hgwdev.cse.ucsc.edu/~benedict/code/Reference_paper/data.zip), expand it and place the data/ directory in the same directory as this README. Then type

    make test

and off you go. There are larger and longer running tests inside of the Makefile if you'd like to justify spending a couple hours doing something else.

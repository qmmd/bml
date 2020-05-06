---
title: BML
---

[![Build Status](https://travis-ci.org/lanl/bml.svg?branch=master)](https://travis-ci.org/lanl/bml)
[![codecov.io](https://codecov.io/github/lanl/bml/coverage.svg?branch=master)](https://codecov.io/github/lanl/bml?branch=master)
[![Coverage Status](https://coveralls.io/repos/github/lanl/bml/badge.svg?branch=master)](https://coveralls.io/github/lanl/bml?branch=master)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/13443/badge.svg)](https://scan.coverity.com/projects/lanl-bml)
[![GitHub issues](https://img.shields.io/github/issues/lanl/bml.svg)](https://github.com/lanl/bml/issues)
[![GitHub pull requests](https://img.shields.io/github/issues-pr/lanl/bml.svg)](https://github.com/lanl/bml/pulls)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/d34f276064ea4421a4d7a3380910cb1c)](https://www.codacy.com/app/nicolasbock/bml?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=lanl/bml&amp;utm_campaign=Badge_Grade)
[![CircleCI](https://circleci.com/gh/lanl/bml/tree/master.svg?style=svg)](https://circleci.com/gh/lanl/bml/tree/master)
![](https://github.com/lanl/bml/workflows/CI/badge.svg)

# Introduction

This website is intended to provide some guidance on how to get and install
the bml library. LA-UR number LA-UR-**17-27373**.

The basic matrix library (bml) is a collection of various matrix data formats
(in dense and sparse) and their associated algorithms for basic matrix
operations. Application programming interfaces (API) are available for both C
and FORTRAN. The current status of this library allows us to use two different
formats for representing matrix data: dense and ELLPACK.

# Mailing List

We are running the following mailing list for discussions on usage and features of the bml library:

<div class="classictemplate template" style="display: block;">
<style type="text/css">
  #groupsio_embed_signup input {border:1px solid #999; -webkit-appearance:none;}
  #groupsio_embed_signup label {display:block; font-size:16px; padding-bottom:10px; font-weight:bold;}
  #groupsio_embed_signup .email {display:block; padding:8px 0; margin:0 4% 10px 0; text-indent:5px; width:58%; min-width:130px;}
  #groupsio_embed_signup {
    background:#fff; clear:left; font:14px Helvetica,Arial,sans-serif;
  }
  #groupsio_embed_signup .button {

      width:25%; margin:0 0 10px 0; min-width:90px;
      background-image: linear-gradient(to bottom,#337ab7 0,#265a88 100%);
      background-repeat: repeat-x;
      border-color: #245580;
      text-shadow: 0 -1px 0 rgba(0,0,0,.2);
      box-shadow: inset 0 1px 0 rgba(255,255,255,.15),0 1px 1px rgba(0,0,0,.075);
      padding: 5px 10px;
      font-size: 12px;
      line-height: 1.5;
      border-radius: 3px;
      color: #fff;
      background-color: #337ab7;
      display: inline-block;
      margin-bottom: 0;
      font-weight: 400;
      text-align: center;
      white-space: nowrap;
      vertical-align: middle;
    }
</style>
<div id="groupsio_embed_signup">
<form action="https://groups.io/g/bml/signup?u=5022678804200651313" method="post" id="groupsio-embedded-subscribe-form" name="groupsio-embedded-subscribe-form" target="_blank">
    <div id="groupsio_embed_signup_scroll">
      <label for="email" id="templateformtitle">Subscribe to our group</label>
      <input type="email" value="" name="email" class="email" id="email" placeholder="email address" required="">

    <div style="position: absolute; left: -5000px;" aria-hidden="true"><input type="text" name="b_5022678804200651313" tabindex="-1" value=""></div>
    <div id="templatearchives"><p><a id="archivelink" href="https://groups.io/g/bml/topics">View Archives</a></p></div>
    <input type="submit" value="Subscribe" name="subscribe" id="groupsio-embedded-subscribe" class="button">
  </div>
</form>
</div>
</div>

# Supported Matrix Formats

The bml library supports the following matrix formats:

* dense
* ELLPACK
* ELLSORT

# Binary Packages

We offer binary packages of the bml library in [RPM
format](http://software.opensuse.org/download.html?project=home%3Anicolasbock%3Aqmmd&package=bml)
thanks to SUSE's OpenBuild Service and for Ubuntu in [DEB
format](https://launchpad.net/~nicolasbock/+archive/ubuntu/qmmd).

# Build Instructions

The bml library is built with CMake. For convenience, we provide a shell
script which goes through the necessary motions and builds the library, runs
the tests, and installs it (in the `install` directory). Simply run:

    $ ./build.sh install

and the library will be built in the `build` directory and installed in the
`install` directory. In case you change any sources and simply want to
rebuild the library, you don't have to run `build.sh` again, but rather

    $ make -C build

The compiled library can be installed by running

    $ make -C build install

The install directory can be modified by running

    $ CMAKE_INSTALL_PREFIX=/some/path ./build.sh install

(which assumes that you are using the bash shell).

To build with GNU compilers, OpenMP, and Intel MKL do the following.

    $ CC=gcc FC=gfortran \
        BLAS_VENDOR=Intel CMAKE_BUILD_TYPE=Release \
        BML_OPENMP=yes CMAKE_INSTALL_PREFIX=/some/path \
        ./build.sh install

To build with MPI, OpenMP, and use Intel MKL do the following.

    $ CC=mpicc FC=mpif90 \
        BLAS_VENDOR=Intel CMAKE_BUILD_TYPE=Release \
        BML_OPENMP=yes BML_MPI=yes CMAKE_INSTALL_PREFIX=/some/path \
        ./build.sh install

## Prerequisites

In order to build the library, the following tools need to be installed:

- `gcc` with Fortran support
- `>=cmake-2.8.8`
- `>=python-2.7`
- `>=OpenMP-3.1` (i.e. `>=gcc-4.7`)

## If the build fails

In case the build fails for some reason, please contact the developers by
opening an issue on GitHub (https://github.com/lanl/bml/issues) and attach the
files

    build/CMakeFiles/CMakeOutput.log
    build/CMakeFiles/CMakeError.log

# Developer Suggested Workflow

Our main development happens on the `master` branch and is continuously
verified for correctness. If you would like to contribute your work to the bml
project, please follow the instructions at the GitHub help page ["About pull
requests"](https://help.github.com/articles/about-pull-requests/). To
summarize:

- Fork the project on github
- Clone that forked repository
- Create a branch in it
- Commit any changes to the branch
- Push the branch to your forked repository
- Go to https://github.com/lanl/bml and click on 'Create Pull Request'

During the review process you might want to update your pull
request. Please add commits or `amend` your existing commits as
necessary. If you amend any commits you need to add the
`--force-with-lease` option to the `git push` command. Please make
sure that your pull request contains only one logical change (see
["Structural split of
change"](https://wiki.openstack.org/wiki/GitCommitMessages#Structural_split_of_changes)
for further details.

## Coding Style

Please indent your C code using

    $ indent -gnu -nut -i4 -bli0 -cli4 -ppi0 -cbi0 -npcs -bfda

You can use the script `indent.sh` to indent all C code.

# Citing

If you find this library useful, we encourage you to cite us. Our project has
a citable DOI:

[![DOI](https://zenodo.org/badge/20454/qmmd/bml.svg)](https://zenodo.org/badge/latestdoi/20454/qmmd/bml)

with the following `bibtex` snipped:

    @misc{bml,
      author       = {Nicolas Bock and
                      Susan Mniszewski and
                      Bálint Aradi and
                      Michael Wall and
                      Christian F. A. Negre and
                      Jamal Mohd-Yusof},
      title        = {qmmd/bml v1.2.3},
      month        = feb,
      year         = 2018,
      doi          = {10.5281/zenodo.841949},
      url          = {https://doi.org/10.5281/zenodo.841949}
    }

# Authors

The core developers of the bml in alphabetical order:

* Christian Negre <cnegre@lanl.gov>
* Nicolas Bock <nicolasbock@gmail.com>
* Susan M. Mniszewski <smm@lanl.gov>

# Contributors

* Jamaludin Mohd-Yusof <jamal@lanl.gov>
* Bálint Aradi <aradi@uni-bremen.de>
* Mike Wall <mewall@lanl.gov>

# License

The bml library is licensed under the BSD 3-clause license.

Copyright 2015. Los Alamos National Security, LLC. This software was
produced under U.S. Government contract DE-AC52-06NA25396 for Los
Alamos National Laboratory (LANL), which is operated by Los Alamos
National Security, LLC for the U.S. Department of Energy. The
U.S. Government has rights to use, reproduce, and distribute this
software. NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY,
LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY
FOR THE USE OF THIS SOFTWARE. If software is modified to produce
derivative works, such modified software should be clearly marked, so
as not to confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with
or without modification, are permitted provided that the following
conditions are met:
- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
- Neither the name of Los Alamos National Security, LLC, Los Alamos
  National Laboratory, LANL, the U.S. Government, nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS
ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# LA-CC

NOTICE OF OSS COPYRIGHT ASSERTION:

LANS has asserted copyright on the software package entitled *Basic
Matrix Library (bml), Version 0.x (C16006)*.

## ABSTRACT

The basic matrix library (bml) is a collection of various matrix data
formats (in dense and sparse) and their associated algorithms for basic
matrix operations.

This code is unclassified and has been assigned LA-CC-**15-093**. Los Alamos
National Laboratory’s Export Control Team made an in-house determination that
this software is controlled under Department of Commerce regulations and the
Export Control Classification Number (ECCN) **EAR99**. The export control
review is attached.

The developers intend to distribute this software package under the OSI
Certified **BSD 3-Clause License**
(http://www.opensource.org/licenses/BSD-3-Clause)

This code was developed using funding from:

- Basic Energy Sciences (LANL2014E8AN) and the Laboratory Directed Research
  and Development Program of Los Alamos National Laboratory. To tests these
  developments we used resources provided by the Los Alamos National
  Laboratory Institutional Computing Program, which is supported by the U.S.
  Department of Energy National Nuclear Security Administration

- Exascale Computing Project (17-SC-20-SC), a collaborative effort of two U.S.
  Department of Energy organizations (Office of Science and the National
  Nuclear Security Administration) responsible for the planning and
  preparation of a capable exascale ecosystem, including software,
  applications, hardware, advanced system engineering, and early testbed
  platforms, in support of the nation’s exascale computing imperative.

Larry Kwei, LAFO Program Manager, has granted his concurrence to asserting
copyright and then distributing the **Basic Matrix Library (bml), Version
0.x** code using an open source software license. See attached memo.

LANS acknowledges that it will comply with the DOE OSS policy as follows:

1. submit form DOE F 241.4 to the Energy Science and Technology Software
   Center (ESTSC),
2. provide the unique URL on the form for ESTSC to distribute, and
3. maintain an OSS Record available for inspection by DOE.

Following is a table briefly summarizes information for this software package:

| CODE NAME                                   | Basic Matrix Library (bml), Version 0.x (C16006) |
| ------------------------------------------- | ------------------------------------------------ |
| Classification Review Number                | **LA-CC-15-093**                                 |
| Export Control Classification Number (ECCN) | **EAR99**                                        |
| B&R Code                                    | **YN0100000**                                    |

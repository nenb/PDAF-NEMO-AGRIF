Quickstart Guide to Building PDAF-NEMO-AGRIF
============================================

 **Step 1** Clone this repository:
``` bash
$ git clone https://github.com/nenb/PDAF-NEMO-AGRIF
```

**Step 2a** Clone the PDAF-OMIv1.2 repository (currently private, hopefully public soon):
``` bash
$ git clone https://github.com/nenb/PDAF_OMIv1.2
```

**Step 2b** Build PDAF-OMIv1.2. See the README file in the PDAF-OMI root directory for build instructions.

**Step 3** Retrieve and build XIOS (https://forge.ipsl.jussieu.fr/ioserver/wiki) on your local machine. This project has been tested with XIOS 2.

**Step 4** In PDAF-NEMO-AGRIF/, modify the appropriate ARCH/ file for your local machine. This will include adding the locations of your newly built PDAF-OMI 
and XIOS executables.

**Step 5** Option 1: In the CONFIG/ subdirectory, build NEMO according to the official documentation 
(https://forge.ipsl.jussieu.fr/nemo/chrome/site/doc/NEMO/guide/html/install.html#compile-and-create-nemo-executable). Make sure to select OPA_SRC, LIM_SRC_2 and
NST_SRC for your build and make sure to select the file with the CPP keys that is relevant to your build (NEMO, NEMO-AGRIF, NEMO-PDAF, NEMO-AGRIF-PDAF are the
current options). Also make sure to include the MY_SRC subdirectory in CONFIG/ in your build.

**Step 5** Option 2: In the CONFIG/ subdirectory, edit the (primitive) shell script build.sh to something appropriate for your local machine. Then run this script
to complete the build.

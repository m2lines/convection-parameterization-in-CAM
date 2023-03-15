# Building CESM on CSD3 Icelake nodes.

- Clone the repo and checkout the release branch you want:
```bash
git clone https://github.com/ESCOMP/CESM
git checkout release-cesm2.1.3
```
- Next, run:
```bash
python3 ./manage_externals checkout_externals
```
which will prompt you to agree to a certificate. You can reject, or accept on temporary or permanent basis. So far, we tried temporary, which fails the first time you run it, yet apparently succeeds the second time (for reasons).
- According to the README, we now have a working version of CESM.
- The README then refers us to [this](http://esmci.github.io/cime/versions/master/html/index.html) documentation to set up a case.

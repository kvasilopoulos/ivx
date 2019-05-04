## Test environments

* local OS MS install, R 3.6.0
* ubuntu 14.04 on travis-ci (devel and release)
* macOS on travis-ci (devel and release)
* r-hub: windows-x86_64-devel, ubuntu-gcc-release, fedora-clang-devel
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* Maintainer: 'Kostas Vasilopoulos <k.vasilopoulo@gmail.com>'

New submission

## General Comments

* Added examples on all exported functions.

* The procedure is called IVX  because the instrumental variable relies directly 
on the regressor, X. Since it is not exactly an acronym it is not explained in the description.

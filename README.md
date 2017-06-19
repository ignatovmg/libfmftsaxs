# libmfftsaxs
### Provides tools for manifold-FFT SAXS calculation.

---

##Requirements
* libmol2
* Check
* Doxygen
* CMake (version > 3.0)

##Installation
After cloning the repository, do
```bash
$ cd libmfftsaxs
$ mkdir build && cd build
$ cmake -DBUILD_TESTS=True -DBUILD_DOCS=True -DCMAKE_BUILD_TYPE=Release ..
$ make
$ make test && make doc
$ make install
```
---
Project is in development.

---
##Coding style
For coding style, refer to libmol2 (especially atom_group.c/h module).  
Except that all interfaced functions should start with **smp\_** instead of **mol\_**:  
Ex: ``void smp_object_action()``

All indentation is with tabs.

Please, follow the style strictly.   



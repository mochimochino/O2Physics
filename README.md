# O2Physics

Documentation for the ALICE Analysis framework is available at:

<https://aliceo2group.github.io/analysis-framework/>

---

## About this branch (Fast Debugging Mode)

This branch has been customized for rapid development and debugging of specific tutorial and analysis tasks, currently focusing on `UDTutorial_06.cxx` and `UDTutorial_07.cxx`.

### Future Plans & Direction
The main goal of this branch is to evolve into a specialized, easy-to-use local environment dedicated to **UD (Ultra-peripheral & Diffractive physics) tasks**. 
By stripping away heavy dependencies and unnecessary modules, this branch will serve as a lightweight building and testing ground, allowing developers to quickly iterate and debug UD-related analysis workflows locally without the overhead of compiling the entire O2Physics repository.

### Current Optimizations
To significantly reduce compilation time (especially when using limited CPU cores like `-j 2`), the `CMakeLists.txt` files have been modified to exclude unnecessary modules and tutorials from the build loop. The following changes were applied:

- **`/CMakeLists.txt`**: Dropped most of the high-level analysis directories (e.g., `PWGCF`, `PWGDQ`, `EventFiltering`, etc.). Only `Common` and `Tutorials` are built.
- **`/Common/CMakeLists.txt`**: Disabled extra framework tools like `Tasks`, `TableProducer`, `Tools`, and `LegacyDataQA` to eliminate missing dependencies on unused modules.
- **`/Tutorials/CMakeLists.txt`**: Skipped building all other tutorial subgroups, isolating only `PWGUD`.
- **`/Tutorials/PWGUD/CMakeLists.txt`**: Enclosed non-target workflows in `if(FALSE)` ... `endif()`, so only required tasks (`udtutorial-06` and `udtutorial-07`) are compiled.

This allows you to quickly re-run `aliBuild build O2Physics` and test changes iteratively.
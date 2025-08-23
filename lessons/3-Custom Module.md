# Custom Module

### Creating a custom module

If you have a program that you want to make available to others, you can create a custom module for it. This involves creating a module file that defines the environment variables and paths needed to run your program.

**Example module file (`myprogram/1.0.lua`):**

```lua
help([[This module loads MyProgram version 1.0]])

prepend_path("PATH", "/path/to/myprogram/bin")
prepend_path("LD_LIBRARY_PATH", "/path/to/myprogram/lib")
```

Place this file in a directory that is part of the `MODULEPATH` environment variable. You can check your `MODULEPATH` with `module use`.

### Using a custom module

Once your custom module is created and placed in the correct location, you can load it like any other module:

```bash
module load myprogram/1.0
```



import subprocess
import os
from hatchling.builders.hooks.plugin.interface import BuildHookInterface

class CustomBuildHook(BuildHookInterface):
    def initialize(self, version, build_data):
        self.compile_dnapars()

    def compile_dnapars(self):
        # Define the paths to your C files
        source_files = ["tribal/dnapars/phylip.c", "tribal/dnapars/seq.c", "tribal/dnapars/dnapars.c"]
        object_files = ["tribal/dnapars/phylip.o", "tribal/dnapars/seq.o", "tribal/dnapars/dnapars.o"]
        executable = "tribal/dnapars/dnapars"

        # Compile the C source files
        for source, obj in zip(source_files, object_files):
            subprocess.run(["cc", "-c", source, "-o", obj], check=True)

        # Link the object files to create the executable
        subprocess.run(["cc", *object_files, "-o", executable, "-lm"], check=True)

        # Make the executable globally accessible (optional)
        os.chmod(executable, 0o755)

"""Generate the code reference pages."""

import os
from pathlib import Path

import mkdocs_gen_files

SRC_DIR = "gypsum_dl"
WRITE_DIR = "api"

for path in sorted(Path(SRC_DIR).rglob("*.py")):  #
    module_path = path.relative_to(SRC_DIR).with_suffix("")  #

    doc_path = path.relative_to(SRC_DIR).with_suffix(".md")  #

    if not os.path.exists(Path(WRITE_DIR)):
        os.mkdir(Path(WRITE_DIR))

    full_doc_path = Path(WRITE_DIR, doc_path)  #

    parts = list(module_path.parts)

    if parts[-1] == "__init__":  #
        parts = parts[:-1]
    elif parts[-1] == "__main__":
        continue

    if len(parts) == 0:
        continue

    with mkdocs_gen_files.open(full_doc_path, "w") as fd:  #
        identifier = ".".join(parts)  #

        print("::: " + identifier, file=fd)  #

    mkdocs_gen_files.set_edit_path(full_doc_path, path)  #

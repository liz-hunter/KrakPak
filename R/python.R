.krakpak_py <- new.env(parent = emptyenv())

.krakpak_python <- function() {

  if (is.null(.krakpak_py$module)) {

    python_path <- system.file(
      "python",
      package = "KrakPak"
    )

    .krakpak_py$module <- reticulate::import_from_path(
      "krakpak",
      path = python_path
    )
  }

  .krakpak_py$module
}

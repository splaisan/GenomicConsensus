This actually failed once because of a missing import, so we might as
well test it.

  $ variantCaller --version
  1.1.0

This will break if the parser setup is messed up.

  $ variantCaller --help >/dev/null

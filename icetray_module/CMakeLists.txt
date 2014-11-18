i3_project(royfit
  PYTHON_DIR python
  PYTHON_DEST icecube/royfit
  DOCS_DIR resources/docs
  )

if (PYTHON_NUMERIC_VERSION LESS 20600)
    colormsg (YELLOW "*** python v2.6 or greater is required... skipping")
else (PYTHON_NUMERIC_VERSION LESS 20600)

i3_test_scripts(resources/test/*.py)

endif(PYTHON_NUMERIC_VERSION LESS 20600)

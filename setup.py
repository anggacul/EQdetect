import setuptools
from distutils.core import setup
from Cython.Build import cythonize
setup(
    ext_modules=cythonize(["./EQdetect/utils/__init__.pyx", "./EQdetect/__init__.pyx", "./EQdetect/core/__init__.pyx","./EQdetect/app.pyx","./EQdetect/utils/report.pyx", "./EQdetect/utils/update_db.pyx", "./EQdetect/utils/vorstat.pyx", "./EQdetect/utils/config.pyx","./EQdetect/core/ipf.pyx", "./EQdetect/core/sourcecal.pyx"]),
    install_requires=[
        'mysql-connector',
    ]
)
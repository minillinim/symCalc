from distutils.core import setup

setup(
    name='symCalc',
    version='0.0.1',
    author='Michael Imelfort',
    author_email='mike@mikeimelfort.com',
    packages=['symcalc', 'symcalc.test'],
    scripts=['bin/symCalc'],
    url='http://pypi.python.org/pypi/symCalc/',
    license='GPLv3',
    description='symCalc',
    long_description=open('README.md').read(),
    install_requires=[],
)


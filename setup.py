from distutils.core import setup

setup(
    name='symCalc',
    version='1',
    author='Michael Imelfort',
    author_email='mike@mikeimelfort.com',
    packages=['symcalc', 'symcalc.test'],
    scripts=['bin/symCalc'],
    url='http://pypi.python.org/pypi/symCalc/',
    license='GPLv3',
    description='Use OTU tables as a reference for calculating read counts for multiple sampled metagenomes',
    long_description=open('README.md').read(),
    install_requires=["biom >= 1.3.1",
                      "ParseM >= 0.0.1",
                      "json >= 2.0"
                      ],
)


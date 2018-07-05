from setuptools import setup

setup(name='helprop',
    version='0.0.1.dev',
    description='Calculate state properties for helium.',
    url='',
    author='Alex Morgan',
    author_email='alexandre.morgan.15@ucl.ac.uk',
    license='BSD 3-clause',
    packages=['helprop'],
    install_requires=[
        'tqdm',
        'attrs'
    ],
    include_package_data=True,
    zip_safe=False)

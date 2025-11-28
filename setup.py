from setuptools import find_packages, setup

setup(
    name="xce",
    version="2.0.1",
    description="GUI tool for XChem",
    author="Tobias Krojer",
    packages=find_packages(),
    package_data={
        'xce': [
            'lib/html_fragments/*.html',
            'web/jscss/css/*',
            'web/jscss/js/*',
            'image/*.png',
        ],
    },
    include_package_data=True,
    install_requires=[],
)

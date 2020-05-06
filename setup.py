import setuptools

setuptools.setup(
    name="bmlscripts",
    version="0.1",
    packages=setuptools.find_packages(),
    entry_points={
        "console_scripts": [
            "convert-template = scripts.convert_template:main"
        ]}
)

from setuptools import setup
def readme():
    with open('README.rst') as f:
        return(f.read())

setup(name='weighted_mutual',
      version = '0.1',
      description = 'calculate weighted mutual information for protein functional correlation based on domain sharing',
      long_description = readme(),
      keywords = 'protein function domain information mutual',
      url = 'http://github.com/algaebrown/weighted_mutual',
      author = 'algaebrown',
      author_email = 'b101102109@tmu.edu.tw',
      license = 'MIT',
      packages = ['weighted_mutual'],
      install_requires = ['pandas', 'numpy'],
      test_suite = 'test',
      zip_safe = False)

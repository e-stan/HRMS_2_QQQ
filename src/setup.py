from setuptools import setup,find_packages

setup(
  name = 'srm_helper',         # How you named your package folder (MyLib)
  packages = ["srm_helper"],   # Chose the same as "name"
  version = '0.0.1',      # Start with a small number and increase it with every change you make
  license='MIT',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = 'Metabolomics software for automated generated of SRMs from HRMS MS/MS data',   # Give a short description about your library
  author = 'Ethan Stancliffe',                   # Type in your name
  author_email = 'estancl1234@gmail.com',      # Type in your E-Mail
  url = 'https://github.com/e-stan/HRMS_2_QQQ',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/e-stan/HRMS_2_QQQ/archive/v0.0.1.tar.gz',    # I explain this later on
  keywords = ['Metabolomics', 'SRM', 'MS/MS',"QqQ"],   # Keywords that define your package best
  install_requires=[            # I get to this in a second
          'numpy',
          'sklearn',
          'pandas',
          'DecoID',
          'matplotlib'
      ],
  classifiers=[
    'Development Status :: 4 - Beta',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Developers',      # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',   # Again, pick a license
    'Programming Language :: Python :: 3.7',      #Specify which pyhton versions that you want to support
  ],

)
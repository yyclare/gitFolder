#!/usr/bin/env python

import os
import ConfigParser


class Config():
    """
        Read element from the configuration file.
    """
    def __init__(self, openedConfigFile, section):
        """
        @param configFile:
        @param section:
        @return:
        """
        self._config = ConfigParser.ConfigParser()
        self._config.readfp(openedConfigFile)
        self._section = section
        self._configFilePath = openedConfigFile.name
        #if defaultDict is not None:
         #   self._defaultDict = defaultDict
        #else:
        #    self._defaultDict = dict()

    def getConfigFile(self):
        return self._configFilePath

    def get(self, option):
        """
            Gets an element 'option' from a given section.
        """
        elem = self._config.get(self._section, option)
        if elem != '':
            return elem
        #else:
        #    return self._defaultDict.get(option, None)


def _test():
    config = Config(open(os.path.normpath('/net/metagenomics/projects/shotgun_pipeline/sources/shotgun_config.cfg')), 'SMP')
    print config.get('inputFastqFileForward1')


if __name__ == "__main__":
    _test()

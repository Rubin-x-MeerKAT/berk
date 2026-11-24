"""

User settings for Berk (e.g., credentials and URLs for optical survey databases).

Settings are resolved in this order:

    1. A value passed explicitly to the function that needs it (keyword argument).
    2. An environment variable (e.g., ``RUBIN_TOKEN``).
    3. The Berk config file - this is ``$BERK_CONFIG`` if that environment variable is set, otherwise
       ``~/.config/berk/config.toml``.
    4. A default value (if there is one).

The config file is in TOML format, with settings grouped into sections, e.g.::

    [crossmatch]
    zcluster_cache = "~/.cache/berk"
    rubin_url = "https://..."
    rubin_token = "..."
    dl_username = "..."
    dl_password = "..."

Keys in the config file are the lower case versions of the corresponding environment variable names.
Since the config file may contain tokens and passwords, it should only be readable by you
(``chmod 600 ~/.config/berk/config.toml``).

Unlike :mod:`berk.startup`, this module has no side effects on import.

"""

import os
import stat
import warnings

try:
    import tomllib
except ModuleNotFoundError:
    try:
        import tomli as tomllib
    except ModuleNotFoundError:
        tomllib = None

DEFAULT_CONFIG_PATH = os.path.join("~", ".config", "berk", "config.toml")

#------------------------------------------------------------------------------------------------------------
def configFilePath():
    """Returns the path to the Berk config file (which may not exist).

    """
    return os.path.expanduser(os.environ.get('BERK_CONFIG', DEFAULT_CONFIG_PATH))

#------------------------------------------------------------------------------------------------------------
def readConfigFile():
    """Reads the Berk config file.

    Returns:
        A dictionary of the config file contents (empty if the config file does not exist).

    """
    path = configFilePath()
    if os.path.exists(path) == False:
        return {}
    if tomllib is None:
        raise Exception("Reading the Berk config file (%s) needs Python >= 3.11, or the tomli package to be installed." % (path))
    with open(path, "rb") as inFile:
        return tomllib.load(inFile)

#------------------------------------------------------------------------------------------------------------
def _checkPermissions(path):
    """Warns if the config file is readable by users other than its owner.

    """
    mode = os.stat(path).st_mode
    if mode & (stat.S_IRWXG | stat.S_IRWXO):
        warnings.warn("Berk config file %s contains secrets but is readable by other users - fix with: chmod 600 %s" % (path, path))

#------------------------------------------------------------------------------------------------------------
def get(envVar, value = None, section = 'crossmatch', default = None, required = True, secret = False):
    """Returns a setting, looking in turn at `value`, the environment variable `envVar`, the config file,
    and finally `default`.

    Args:
        envVar (:obj:`str`): Name of the environment variable (e.g., 'RUBIN_TOKEN'). The corresponding
            key in the config file is the lower case version of this (e.g., 'rubin_token').
        value (optional): A value given explicitly (e.g., by a keyword argument). If not None, this is
            returned as is.
        section (:obj:`str`, optional): Section of the config file in which to look for the setting.
        default (optional): Value to return if the setting is not found anywhere else.
        required (:obj:`bool`, optional): If True, and the setting is not found (and `default` is None),
            raise an exception.
        secret (:obj:`bool`, optional): If True, warn if the value is read from a config file that is
            readable by other users.

    Returns:
        The value of the setting (None if not found, and `required` is False).

    """

    if value is not None:
        return value

    if envVar in os.environ.keys():
        return os.environ[envVar]

    key = envVar.lower()
    config = readConfigFile()
    if section in config.keys() and key in config[section].keys():
        if secret == True:
            _checkPermissions(configFilePath())
        return config[section][key]

    if default is not None:
        return default

    if required == True:
        raise Exception("Setting '%s' not found - pass it as a keyword argument, set the %s environment variable, or add '%s = ...' to the [%s] section of the Berk config file (%s)."
                        % (key, envVar, key, section, configFilePath()))

    return None

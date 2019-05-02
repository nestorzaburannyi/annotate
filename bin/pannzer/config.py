from __future__ import print_function
import sys
from os.path import isfile
if list(sys.version_info)[0] == 2:
 from ConfigParser import ConfigParser
else:
 from configparser import ConfigParser

class BugsyConfig(object) :
    def __init__(self, fname) :
        self.fname = fname
        self.config = ConfigParser.ConfigParser()

        if isfile(self.fname) :
            self.config.read(self.fname)

    def get(self, section, option) :
        try :
            return self.config.getint(section, option)
        except ValueError :
            pass

        try :
            return self.config.getfloat(section, option)
        except ValueError :
            pass

        try :
            return self.config.getboolean(section, option)
        except ValueError :
            pass

        return self.config.get(section, option)

    def set(self, section, option, value) :
        assert None not in (section, option, value)

        try :
            self.config.set(section, option, str(value))

        except ConfigParser.NoSectionError:
            self.config.add_section(section)
            self.config.set(section, option, str(value))

        print("set %s::%s to %s" % (section, option, value))

    def flush(self, fname=None) :
        if not fname :
            fname = self.fname

        with open(fname, 'w') as f :
            self.config.write(f)

    def echo(self, section=None, option=None) :
        if section == None :
            for s in self.config.sections() :
                print("[%s]" % s)
                for k,v in self.config.items(s) :
                    print("  %s = %s" % (k,v))

        elif option == None :
            if not self.config.has_section(section) :
                print("Error: section '%s' does not exist" % section, file=stderr)
                exit(1)

            for k,v in self.config.items(section) :
                print("%s = %s" % (k,v))

        else :
            if not self.config.has_section(section) :
                print("Error: section '%s' does not exist" % section, file=stderr)
                exit(1)

            if not self.config.has_option(section, option) :
                print("Error: option '%s' does not exist in section '%s'" % (option, section), file=stderr)
                exit(1)

            print(self.config.get(section, option))

def usage() :
    print("Usage:\t%s set <fname> <section> <option> <value>\n\t%s get <fname> [<section> [<option>]]" % (argv[0], argv[0]), file=stderr)

def main() :
    if len(argv) < 3 :
        usage()
        return 1
    
    command  = argv[1]
    fname    = argv[2]

    if command not in ('set', 'get') :
        usage()
        return 0

    if command == 'set' :
        if len(argv) < 6 :
            usage()
            return 1

    section = option = value = None

    try :
        section = argv[3]
        option  = argv[4]
        value   = argv[5]

    except IndexError :
        pass

    bc = BugsyConfig(fname)

    if command == 'set' :
        bc.set(section, option, value)
        bc.flush()

    elif command == 'get' :
        if not isfile(fname) :
            print("Error: '%s' not found" % fname, file=stderr)
            return 1

        bc.echo(section, option)

    return 0


if __name__ == '__main__' :
    try :
        exit(main())

    except KeyboardInterrupt :
        exit(1)

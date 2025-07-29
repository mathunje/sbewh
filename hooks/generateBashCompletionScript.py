#!/usr/bin/env python3


import argparse
import subprocess
import sys


def main(sbewhPath, optionsPath=None):
    if not optionsPath is None:
        options = [ "-" + s[:-1] for s in open(optionsPath, "r").readlines() if len(s) > 0 ]
    else:
        result = subprocess.run([sbewhPath + " -par"], shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Call of sbewh ({sbewhPath -par}) failed")
        options = [ "-" + s for s in result.stdout.split("\n") if len(s) > 0 ]

    groups = {}
    for opt in options:
        g = opt.partition(".")[0]
        if g in groups:
            groups[g].append(opt)
        else:
            groups[g] = [opt]
    additionalKeys = ["-h", "-help", '-par']
    additionalKeysArr = "("  + " ".join([f"\"{s}\"" for s in additionalKeys]) + ")"
    allKeys = list(groups.keys()) + additionalKeys
    mainGroupCompReply = f"COMPREPLY=( $(compgen -W \"{" ".join(allKeys)}\" -- $cword) )"
    subMatch = "      case ${group} in\n"
    for group, opts in groups.items():
        subMatch += f"        {group})\n"
        subMatch += f"            COMPREPLY=( $(compgen -W \"{" ".join(opts)}\" -- $cword) );;\n"
    subMatch += "  esac\n"


    script="""#/usr/bin/env bash

_sbewh_completions()
{
    local additionalKeys=""" + additionalKeysArr + """
    cword="${COMP_WORDS[$COMP_CWORD]}"
    helpMode=0
    if [ $COMP_CWORD -gt 0 ]; then
        prev=${COMP_WORDS[$COMP_CWORD-1]}
        if ( [ "${prev}" = "-help" ] || [ "${prev}" = "-h" ]) ; then
            helpMode=1
            if  [ "${cword:0:1}" != "-" ] ; then
                cword="-${cword}"
            fi
        fi
    fi
    if [ "${cword:0:1}" != "-" ] ; then
        if [ ${helpMode} -eq 0 ] ; then
            COMPREPLY=( $(compgen -f -- "$cword"; compgen -d -S / -- "$cword") )
        fi
        return
    fi
    group=$(echo ${cword} | awk -F '.' '{print $1}')
    if [ ${#group} -eq ${#cword} ] ; then
        """ + mainGroupCompReply + """
    fi
    if [ ${#COMPREPLY[@]} -eq 1 ] ; then
        if [[ ${additionalKeys[@]} =~ $cword ]] ; then
            return
        fi
        if [ ${helpMode} -eq 0 ] ; then
            cword="${COMPREPLY[0]}."
        fi
    fi
    group=$(echo ${cword} | awk -F '.' '{print $1}')
""" + subMatch + """
    for key in "${additionalKeys[@]}" ; do
        if [[ ${COMP_WORDS[@]} =~ "${key}" ]] ; then
            return
        fi
    done
    if [ ${#COMPREPLY[@]} -eq 1 ] && [ ${helpMode} -eq 0 ] ; then
        COMPREPLY=( "${COMPREPLY[0]}=")
    fi
}

complete -F _sbewh_completions -o nospace ./sbewh"""
    print(script)


def createParser():
    parser = argparse.ArgumentParser(
                        description="""prints a bash script for autocompletion of sbweh excutable (source it)
                                       You may want to adapt the completion path (last line of the script)""",
                        epilog="from MT")
    parser.add_argument('-e', '--sbewhPath', help='path of sbewh execuable', default="../build/sbewh")
    parser.add_argument('-o', '--optionsPath', help='file name with options separated line by line')
    return parser


if __name__ == "__main__":
    parser = createParser()
    args = parser.parse_args()
    main(**vars(args))

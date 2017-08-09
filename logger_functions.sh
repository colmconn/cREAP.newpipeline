#!/bin/bash

## \e doesn't work on macOS terminal so replace it with \033
## see
## http://apple.stackexchange.com/questions/74777/echo-color-coding-stopped-working-in-mountain-lion

green='\033[0;37;42m'
yellow='\033[0;37;43m'
red='\033[0;37;41m'
endColor='\033[0m'

prefix_tag="***"

#enable echo

function message {
    echo "${prefix_tag} ${message}"
}

function info_message {
    local message="$1"
    if [[ -t 1 ]] ; then
	echo -ne "${green}${prefix_tag}${endColor} ${message}"
    else
	echo -n "${prefix_tag} ${message}"
    fi
}


function warn_message {
    local message="$1"
    if [[ -t 1 ]] ; then
	echo -ne "${yellow}${prefix_tag}${endColor} ${message}"
    else
	echo -n "${prefix_tag} ${message}"
    fi
}

function error_message {
    local message="$1"
    if [[ -t 1 ]] ; then
	echo -ne "${red}${prefix_tag}${endColor} ${message}"
    else
	echo -n "${prefix_tag} ${message}"
    fi
}

function info_message_ln {
    local message="$1"
    if [[ -t 1 ]] ; then
	echo -e "${green}${prefix_tag}${endColor} ${message}"
    else
	echo "${prefix_tag} ${message}"
    fi
}


function warn_message_ln {
    local message="$1"
    if [[ -t 1 ]] ; then
	echo -e "${yellow}${prefix_tag}${endColor} ${message}"
    else
	echo "${prefix_tag} ${message}"
    fi
}

function error_message_ln {
    local message="$1"
    if [[ -t 1 ]] ; then
	echo -e "${red}${prefix_tag}${endColor} ${message}"
    else
	echo "${prefix_tag} ${message}"
    fi
}

# info_message  "This is a no new line test info message"
# warn_message  "This is a no new line test warning message"
# error_message "This is a no new line test error message"
# echo
# info_message_ln  "This is a test info message with new line"
# warn_message_ln  "This is a test warning message with new line"
# error_message_ln "This is a test error message with new line"

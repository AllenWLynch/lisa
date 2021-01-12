import lisa
from lisa import FromRegions, FromGenes, FromCoverage
import argparse
from lisa.cli.cli import oneshot_parser, multi_parser, regions_parser, coverage_parser
import argparse
import os
from argparse import HelpFormatter
from lisa.cli.cli import RstFormatter
from lisa.core.utils import Log

def new_format_help(self):

    formatter = RstFormatter(prog = self.prog)

    # description
    formatter.add_text(self.description)

    # usage
    formatter.add_usage(self.usage, self._actions,
                        self._mutually_exclusive_groups)

    # positionals, optionals and user-defined groups
    for action_group in self._action_groups:
        formatter.start_section(action_group.title)
        formatter.add_text(action_group.description)
        formatter.add_arguments(action_group._group_actions)
        formatter.end_section()

    # epilog
    formatter.add_text(self.epilog)

    # determine help from format above
    return formatter.format_help()

def make_cli_page(parsers):

    docs = lisa.cli.cli.__doc__

    docs += '\n\n'.join([new_format_help(m) for m in parsers]).replace('make_docs.py', 'lisa')

    return docs


def make_api_page(modules):

    docs = lisa.__doc__

    docs += '\n\n'.join([m.get_docs() for m in modules])

    return docs

if __name__ == "__main__":
    
    api_docs = make_api_page([FromGenes, FromRegions, FromCoverage])

    with open('python_api.rst', 'w') as f:
        print(api_docs, file = f)

    cli_docs = make_cli_page([oneshot_parser, multi_parser, regions_parser, coverage_parser])

    with open('cli.rst', 'w') as f:
        print(cli_docs, file = f)
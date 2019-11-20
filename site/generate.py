#!/usr/bin/env python
"""
This is the website source code for the paper "Epigenetic changes in aging human monocytes".
See: https://artyomovlab.wustl.edu/aging/

NOTE: python3 required
> source activate py3.5

author oleg.shpynov@jetbrains.com
"""
import datetime
import os
import re
import shutil

from collections import namedtuple

# Tools versions
SPAN_BUILD = '0.11.0.4882'
SPAN_DATE = 'May 17, 2019'

JBR_BUILD = '1.0.beta.4882'
JBR_DATE = 'May 17, 2019'

DIR = os.path.dirname(__file__)
OUT_FOLDER = os.path.join(DIR, 'out')

DistrDescriptor = namedtuple('DistrDescriptor', ['title', 'suffix', 'folder'])


def generate_jbr_data():
    template = 'jetbrains_research_jbr.html'
    with open(os.path.join(DIR, template), 'r') as f:
        template_html = f.read()

    def create_tr(dd, build):
        code_base = "https://download.jetbrains.com/biolabs/jbr_browser"
        fname = "jbr-{}{}".format(build, dd.suffix)
        url = "{}/{}/{}".format(code_base, dd.folder, fname)
        return """
        <tr>
            <td> <a href="{}">{}</a></td>
            <td>{}</td>
        </tr>
        """.format(url, fname, dd.title)

    with open(os.path.join(OUT_FOLDER, template), 'w') as f:
        descrs = [
            DistrDescriptor("Windows 64-bit ZIP archive (includes bundled 64-bit Java Runtime)",
                            "_x64.zip", "win"),
            DistrDescriptor("Windows 32-bit ZIP archive (includes bundled 32-bit Java Runtime)",
                            "_x86.zip", "win"),
            DistrDescriptor("Mac installer (includes bundled 64-bit Java Runtime)",
                            ".dmg", "mac"),
            DistrDescriptor("Linux archive (includes bundled 64-bit Java Runtime)",
                            ".tar.gz", "linux"),
        ]

        content = template_html. \
            replace('@TABLE@', '\n'.join([create_tr(d, JBR_BUILD) for d in descrs])) \
            .replace('@BUILD@', JBR_BUILD) \
            .replace('@DATE@', JBR_DATE)

        seen_file_names = set()
        for dd in descrs:
            fname = '@FILENAME-{}@'.format(dd.folder)
            if fname not in seen_file_names:
                seen_file_names.add(fname)
                content = content.replace(fname, "{}{}".format(JBR_BUILD, dd.suffix))

        f.write(content)


def generate_span_data():
    template = 'jetbrains_research_span.html'
    with open(os.path.join(DIR, template), 'r') as f:
        template_html = f.read()

    def create_tr(build):
        return """
        <tr>
            <td> <a href="https://download.jetbrains.com/biolabs/span/{0}">{0}</a></td>
            <td> Multi-platform JAR package </td>
        </tr>
        """.format("span-{0}.jar".format(build))

    with open(os.path.join(OUT_FOLDER, template), 'w') as f:
        f.write(template_html.
                replace('@TABLE@', create_tr(SPAN_BUILD)).
                replace('@BUILD@', SPAN_BUILD).
                replace('@DATE@', SPAN_DATE)
                )


def _cli():
    if os.path.exists(OUT_FOLDER):
        shutil.rmtree(OUT_FOLDER)
    os.mkdir(OUT_FOLDER)

    print('Creating JetBrains Research JBR page')
    generate_jbr_data()

    print('Creating JetBrains Research SPAN page')
    generate_span_data()

    print('Done')


if __name__ == "__main__":
    _cli()

"""Tiered argparse help: hide advanced flags from default ``--help``.

Two-tier CLI surface for ``speconsense`` and ``speconsense-summarize``:
the **default** ``--help`` shows the user-facing flags (Common / Tuning),
and the **expanded** ``--help-advanced`` adds pre-tuned MCL internals,
MAD knobs, and integral-phase disable flags. Without this split the
help wall makes the actually-tunable surface hard to find.

Usage from a CLI module::

    parser = argparse.ArgumentParser(...)
    install_advanced_help(parser)

    advanced = parser.add_argument_group(
        "Advanced (pre-tuned — rarely needed)",
        "Hidden from --help. View with --help-advanced.",
    )
    add_advanced_argument(advanced, "--inflation", type=float, default=4.0,
                          help="MCL inflation parameter (default: 4.0)")
    # ... more advanced flags ...

The Advanced group's flags are stored with their real help string but
present ``argparse.SUPPRESS`` to the default help formatter. The
``--help-advanced`` action restores the strings, prints the full help,
and exits.
"""

from __future__ import annotations

import argparse
from typing import Any


_ADVANCED_HELP_ATTR = "_advanced_help"


class _TieredHelpFormatter(argparse.HelpFormatter):
    """Suppress argument groups whose actions are all hidden.

    argparse normally renders a group's title + description even when
    every action under it has ``help=SUPPRESS``. We override
    ``add_arguments`` so a group with no visible actions ends up with no
    visible heading either.
    """

    def add_arguments(self, actions):
        visible = [a for a in actions if a.help is not argparse.SUPPRESS]
        if not visible and actions:
            # All actions are suppressed — drop the section heading too.
            # The current section was opened by start_section() in
            # _format_action_group; emptying its items list does the job.
            self._current_section.items = []
            return
        super().add_arguments(visible)


class _HelpAdvancedAction(argparse.Action):
    """``--help-advanced``: render every flag including the Advanced tier."""

    def __init__(self, option_strings, dest, default=argparse.SUPPRESS, help=None):
        super().__init__(
            option_strings=option_strings,
            dest=dest,
            nargs=0,
            default=default,
            help=help,
        )

    def __call__(self, parser, namespace, values, option_string=None):
        for action in parser._actions:
            stashed = getattr(action, _ADVANCED_HELP_ATTR, None)
            if stashed is not None:
                action.help = stashed
        # Re-render. The trailing pointer line is rendered by epilog,
        # which we suppress here so it isn't shown in the expanded view.
        original_epilog = parser.epilog
        parser.epilog = None
        try:
            parser.print_help()
        finally:
            parser.epilog = original_epilog
        parser.exit()


def install_advanced_help(parser: argparse.ArgumentParser) -> None:
    """Wire up ``--help-advanced`` and the trailing pointer line.

    Adds a ``--help-advanced`` action to the parser, installs the
    ``_TieredHelpFormatter`` (which hides argument groups whose every
    action is suppressed), and sets an epilog that points users to
    ``--help-advanced``. Call once per parser, before adding the
    Advanced group's arguments.
    """
    parser.formatter_class = _TieredHelpFormatter
    parser.add_argument(
        "--help-advanced",
        action=_HelpAdvancedAction,
        help="Show full help including pre-tuned / internal flags",
    )
    pointer = "For pre-tuned/internal flags, use --help-advanced."
    if parser.epilog:
        parser.epilog = f"{parser.epilog}\n\n{pointer}"
    else:
        parser.epilog = pointer


def add_advanced_argument(
    group: argparse._ArgumentGroup,
    *args: Any,
    help: str = None,
    **kwargs: Any,
) -> argparse.Action:
    """Add an argument that is hidden from default ``--help``.

    The flag is fully functional — parsing and downstream behavior are
    unchanged. Only its rendered help line is suppressed unless the user
    passes ``--help-advanced``. The real ``help`` string is stashed on
    the action for ``_HelpAdvancedAction`` to restore.
    """
    action = group.add_argument(*args, help=argparse.SUPPRESS, **kwargs)
    setattr(action, _ADVANCED_HELP_ATTR, help)
    return action

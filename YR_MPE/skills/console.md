# YR-MPE Console

## Purpose

Open a shell preconfigured with the bundled phylogenetic tools, or create a
PhyloBOT project directory.

## GUI features

- **Startup a new shell** opens a terminal with the plugin's tool directories
  added to `PATH`.
- **Create new PhyloBOT project** creates a project scaffold for command-line
  work.
- The dialog lists the tool paths and commands that will be available.

There is no analysis command line of its own; the console is a launcher.

## Workflow

1. Open the Console plugin.
2. Review the tool PATH list.
3. Start a shell, or create a PhyloBOT project.
4. Run bundled tools by name in the opened shell.

## Pitfalls

- The shell inherits the environment at launch; tools added later require a
  new shell.
- The console does not create workspace items or provenance records; use the
  plugins for tracked analyses.
- Closing the dialog does not necessarily close an already-opened terminal.
- Tool paths come from `config.json`; a missing path is shown but cannot be
  fixed from this dialog.

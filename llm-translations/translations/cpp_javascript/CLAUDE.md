# Project Instructions

## Scope

This project's working directory is fully self-contained, with three
narrow, explicit exceptions listed below. Outside of those exceptions,
never read, list, search, write, or otherwise reference anything outside
this folder — no other parent directories, sibling project folders, or the
home directory. Do not use `ls`, `find`, `cat`, `grep`, or any other shell
command to explore anything above this folder beyond what is explicitly
allowed here.

### Allowed exceptions

1. **Nix-provided compilers** (pre-resolved environment, `../../nix-env.sh`)
   - You may read `../../nix-env.sh`, a pre-resolved snapshot of the Nix
     dev environment (PATH, env vars, etc.) generated outside this session
     with `nix print-dev-env`.
   - To build or run code that needs a compiler not present on this
     machine (e.g. Rust), source it before the command:
     `source ../../nix-env.sh && cargo build --release`
     `source ../../nix-env.sh && ./target/release/<binary>`
   - Never call `nix`, `nix develop`, or any other nix subcommand directly.
     The sandbox cannot create the Unix domain socket the Nix daemon
     requires, so these commands will always fail here. Use the
     pre-resolved `nix-env.sh` instead.
   - If `nix-env.sh` seems outdated or missing a tool, tell the user it
     may need to be regenerated outside this session — do not attempt to
     invoke Nix yourself.
   - Do not read, list, or explore anything else under `../../`.

2. **Distance matrix** (one directory up, `../distance_matrix`)
   - This is written by `../mlp-instances/load` as a side effect of
     loading an instance (see above) — this is the only sanctioned way to
     modify it.
   - You may read files under `../distance_matrix` (e.g. to check the
     result of a load).
   - Never edit or write to `../distance_matrix` directly yourself.

No other files or folders under `../` or `../../` are accessible. If
something else seems missing, ask the user instead of looking for it
elsewhere.

## Execution style

- Do exactly what is asked. Do not propose alternative approaches,
  refactors, architecture changes, or "improvements" that were not
  requested.
- Do not suggest additional features, optimizations, or next steps unless
  explicitly asked for them.
- If a request is ambiguous, ask a clarifying question instead of guessing
  or filling gaps with your own design decisions.
- Keep changes minimal and scoped strictly to the task described.

## Communication

- Report what you did factually and concisely.
- Do not add unsolicited opinions, caveats about "best practices," or
  recommendations beyond what was asked.

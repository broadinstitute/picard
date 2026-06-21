# Contributing to Picard

This document describes how to contribute to the Picard codebase. 
It applies to code, tests, documentation, and build or CI changes.

## Scope of contributions

We welcome contributions that improve Picard, including:

- Bug fixes
- Performance improvements
- New tools or tool extensions
- Test coverage improvements
- Documentation updates
- Build, CI, or release-related improvements

Large or behavior-changing contributions should be discussed in an issue before implementation.

## Issues and discussion

- Use GitHub Issues to report bugs or propose changes.
- Include enough detail to reproduce bugs (command line, inputs, versions).
- General usage questions should be asked on https://bioinformatics.stackexchange.com (with tag 'picard').

## Development workflow

1. Fork the repository.
2. Create a topic branch from `master`.
3. Make focused, logically separated commits.
4. Ensure tests pass locally.
5. Open a pull request against `master`.

Pull requests should be limited in scope. 
Large refactors or design changes may be rejected without prior discussion.

## Code standards

- Picard is written in Java.
- Follow existing code style and patterns in the codebase.
- Prefer clarity over cleverness.
- Public APIs and command-line behavior should remain backward-compatible unless explicitly discussed.

## Tests

- New features and bug fixes **must** include tests.
- Tests should fail before the fix and pass after.
- Avoid introducing unnecessary test dependencies or large test data.

## Documentation

- Public-facing changes should include documentation updates.
- Command-line behavior changes should be reflected in tool documentation.
- Keep documentation concise and factual.

## Build and verification

Before submitting a pull request, run:

```bash
./gradlew test
```
Pull requests that fail CI will not be reviewed until fixed.

## Review process

All changes require review before merging.
Maintainers may request changes for correctness, clarity, scope, or maintainability.
Approval does not guarantee immediate merge.

## Licensing

By contributing to Picard, you agree that your contributions will be licensed under the project’s MIT License. 
See LICENSE.txt for details.

## Attribution

GitHub automatically tracks code contributions. 
Additional attribution requests should be discussed in an issue.
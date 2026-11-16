# tests/test_g_relative_root.py
# A relative -C path must not be resolved twice once a submodule is populated.
import subprocess

import pytest

from tests.utils_for_tests import normalize_whitespace


# Only run this test over the single 'test_submodule' fixture
@pytest.fixture(name="shared_repos")
def fixture_test_submodule_only(get_all_repos):
    return next(
        repo for repo in get_all_repos if repo["submodule_name"] == "test_submodule"
    )


def write_shared_gitmodules(repo_dir, shared_repos):
    """Write the shared .gitmodules, which has 'fxtag' set, into repo_dir."""
    (repo_dir / ".gitmodules").write_text(shared_repos["gitmodules_content"])


# Override the shared `test_repo` fixture to have the proper .gitmodules, with fxtag set:
# without it, update and status exit 0 but don't do any work, just report that fxtag isn't
# present.
@pytest.fixture(name="test_repo")
def fixture_test_repo_with_fxtag(test_repo, shared_repos):
    write_shared_gitmodules(test_repo, shared_repos)
    return test_repo


# These tell us that the update actually did work
@pytest.fixture(name="printed_lines")
def fixture_printed_lines(shared_repos):
    return {
        "update": f"{shared_repos['submodule_name']} up to date.",
        "status": shared_repos["status2"],
    }


@pytest.mark.parametrize("action", ["update", "status"])
def test_relative_C(git_fleximod, test_repo, shared_repos, printed_lines, action):
    submodule_path = shared_repos["subrepo_path"]
    # main test: initial population of the submodule, also with a relative -C
    result = git_fleximod(test_repo, "-C . update")
    assert result.returncode == 0
    assert (test_repo / submodule_path / ".git").exists()

    # secondary test: make sure -C properly tries to re-update the same module
    result = git_fleximod(test_repo, f"-C . {action}")
    assert result.returncode == 0
    # A doubly-resolved relative path used to raise this error.
    assert "NoSuchPathError" not in result.stderr
    assert printed_lines[action] in normalize_whitespace(result.stdout)
    assert (test_repo / submodule_path / ".git").exists()


def test_relative_C_in_worktree(
    git_fleximod, test_repo, shared_repos, printed_lines, tmp_path
):
    submodule_path = shared_repos["subrepo_path"]
    # main test: initial population of the submodule in the main checkout aka `test_repo`
    result = git_fleximod(test_repo, "update")
    assert result.returncode == 0

    # Create worktree for that repo
    worktree = tmp_path / "wt"
    subprocess.run(
        ["git", "worktree", "add", str(worktree)],
        cwd=test_repo,
        check=True,
        capture_output=True,
    )
    # bring the shared .gitmodules into the worktree
    write_shared_gitmodules(worktree, shared_repos)

    # main test: initial update works with relative -C
    result = git_fleximod(worktree, "-C . update")
    assert result.returncode == 0
    assert "NoSuchPathError" not in result.stderr
    submodule_git = worktree / submodule_path / ".git"
    assert submodule_git.exists()
    # st_mtime_ns is the file's modification time in integer nanoseconds, exact enough to compare.
    populated_at = submodule_git.stat().st_mtime_ns

    # Secondary test: re-updating correctly re-updates the module (which is a no-op)
    result = git_fleximod(worktree, "-C . update")
    assert result.returncode == 0
    # Make sure the post-checkout block ran
    assert printed_lines["update"] in normalize_whitespace(result.stdout)
    # Confirm the update was a no-op: populating writes the submodule's .git, so an
    # unchanged mtime means nothing was re-cloned.
    assert submodule_git.stat().st_mtime_ns == populated_at


def test_absolute_C_from_other_cwd(
    git_fleximod, test_repo, shared_repos, printed_lines, tmp_path
):
    submodule_path = shared_repos["subrepo_path"]
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    # main test: initial population of the submodule, with an absolute -C from an arbitrary cwd
    result = git_fleximod(elsewhere, f"-C {test_repo} update")
    assert result.returncode == 0
    assert (test_repo / submodule_path / ".git").exists()

    # Secondary test: re-update that directory, still from arbitrary cwd
    result = git_fleximod(elsewhere, f"-C {test_repo} update")
    assert result.returncode == 0
    # Only the post-checkout block prints "up to date"; it used to be skipped when the cwd
    # is not the repo.
    assert printed_lines["update"] in normalize_whitespace(result.stdout)

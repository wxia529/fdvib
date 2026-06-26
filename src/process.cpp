#include "fdvib.hpp"

#include <cerrno>
#include <cstdio>
#include <stdexcept>
#include <sys/wait.h>
#include <unistd.h>

namespace fdvib {

std::string shell_quote(const std::string &s) {
    std::string q = "'";
    for (char c : s) q += (c == '\'') ? "'\\''" : std::string(1, c);
    return q + "'";
}

int shell_run(const std::string &cmd, const fs::path &cwd, const fs::path &stdout_path) {
    const pid_t pid = fork();
    if (pid < 0) throw std::runtime_error("fork failed");
    if (pid == 0) {
        if (chdir(cwd.c_str()) != 0) _exit(126);
        FILE *f = std::fopen(stdout_path.c_str(), "w");
        if (!f) _exit(126);
        dup2(fileno(f), STDOUT_FILENO);
        dup2(fileno(f), STDERR_FILENO);
        std::fclose(f);
        execl("/bin/sh", "sh", "-c", cmd.c_str(), static_cast<char *>(nullptr));
        _exit(127);
    }
    int status = 0;
    while (waitpid(pid, &status, 0) < 0 && errno == EINTR) {}
    if (WIFEXITED(status)) return WEXITSTATUS(status);
    return 128 + (WIFSIGNALED(status) ? WTERMSIG(status) : 0);
}

} // namespace fdvib

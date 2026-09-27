// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

package main

import (
	"fmt"
	"os"
	"path/filepath"

	"github.com/wxia529/fdvib/internal/diagnostics"
	"github.com/wxia529/fdvib/internal/export"
	initcmd "github.com/wxia529/fdvib/internal/init"
	"github.com/wxia529/fdvib/internal/run"
	"github.com/wxia529/fdvib/internal/settings"
)

// version is injected at build time with
// -ldflags "-X main.version=<ver>"; it mirrors the C++ FDVIB_VERSION.
var version = "unknown"

func usage(out *os.File) {
	fmt.Fprintln(out, "Usage:")
	fmt.Fprintln(out, "  fdvib init {local|gas}")
	fmt.Fprintln(out, "  fdvib -inp fdvib.in")
	fmt.Fprintln(out, "  fdvib modes RESULTS_DIR")
	fmt.Fprintln(out, "  fdvib shm RESULTS_DIR")
	fmt.Fprintln(out, "  fdvib thermo RESULTS_DIR -inp thermo.in")
	fmt.Fprintln(out, "  fdvib --help")
	fmt.Fprintln(out, "  fdvib --version")
	fmt.Fprintln(out, "Options: --debug [--debug-log FILE] (before or after command)")
}

func absolute(p string) string {
	abs, err := filepath.Abs(p)
	if err != nil {
		return p
	}
	return abs
}

func main() { os.Exit(execute(os.Args)) }

func debugArgs(args []string) ([]string, bool, string, error) {
	clean := []string{args[0]}
	enabled, pathSeen := false, false
	path := ""
	for i := 1; i < len(args); i++ {
		switch args[i] {
		case "--debug":
			if enabled {
				return nil, false, "", fmt.Errorf("duplicate --debug option")
			}
			enabled = true
		case "--debug-log":
			if pathSeen {
				return nil, false, "", fmt.Errorf("duplicate --debug-log option")
			}
			pathSeen = true
			if i+1 >= len(args) || len(args[i+1]) == 0 || args[i+1][0] == '-' {
				return nil, false, "", fmt.Errorf("--debug-log requires a file path")
			}
			i++
			path = args[i]
		default:
			clean = append(clean, args[i])
		}
	}
	if pathSeen && !enabled {
		return nil, false, "", fmt.Errorf("--debug-log requires --debug")
	}
	return clean, enabled, path, nil
}

func execute(raw []string) int {
	args, debug, logPath, parseErr := debugArgs(raw)
	if parseErr != nil {
		fmt.Fprintf(os.Stderr, "fdvib: error: %s\n", parseErr)
		return 2
	}
	if len(args) == 2 && args[1] == "--help" {
		usage(os.Stdout)
		return 0
	}
	if len(args) == 2 && args[1] == "--version" {
		fmt.Printf("fdvib %s\n", version)
		return 0
	}
	if len(args) < 3 {
		usage(os.Stderr)
		return 2
	}
	cmd := args[1]
	valid := (cmd == "init" && len(args) == 3) || (cmd == "-inp" && len(args) == 3) || ((cmd == "modes" || cmd == "shm") && len(args) == 3) || (cmd == "thermo" && len(args) == 5 && args[3] == "-inp")
	if !valid {
		usage(os.Stderr)
		return 2
	}
	var err error
	var log *diagnostics.Logger
	if debug {
		log, err = diagnostics.Open(logPath, version, os.Stderr)
		if err != nil {
			fmt.Fprintf(os.Stderr, "fdvib: error: %s\n", err)
			return 1
		}
	}
	defer func() { log.Close(err) }()
	log.Event("command", diagnostics.Fields{"command": cmd})
	switch {
	case cmd == "init" && len(args) == 3:
		cwd, getwdErr := os.Getwd()
		if getwdErr != nil {
			err = getwdErr
			break
		}
		err = log.Do("init", nil, func() error { return initcmd.InitializeInput(args[2], cwd) })
	case cmd == "-inp" && len(args) == 3:
		var s *settings.Settings
		err = log.Do("config.load", diagnostics.Fields{"path": args[2]}, func() error {
			var e error
			s, e = settings.From(args[2], "")
			return e
		})
		if err == nil {
			err = run.CalculateWithDiagnostics(s, log)
		}
	case cmd == "modes" && len(args) == 3:
		err = log.Do("export.modes", nil, func() error { return export.Modes(absolute(args[2])) })
	case cmd == "shm" && len(args) == 3:
		err = log.Do("export.shm", nil, func() error { return export.Shm(absolute(args[2])) })
	case cmd == "thermo" && len(args) == 5 && args[3] == "-inp":
		err = log.Do("export.thermo", nil, func() error { return export.Thermo(absolute(args[2]), absolute(args[4])) })
	default:
		usage(os.Stderr)
		return 2
	}
	if err != nil {
		fmt.Fprintf(os.Stderr, "fdvib: error: %s\n", err)
		return 1
	}
	return 0
}

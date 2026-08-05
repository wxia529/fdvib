// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

package main

import (
	"fmt"
	"os"
	"path/filepath"

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
}

func absolute(p string) string {
	abs, err := filepath.Abs(p)
	if err != nil {
		return p
	}
	return abs
}

func main() {
	args := os.Args
	if len(args) == 2 && args[1] == "--help" {
		usage(os.Stdout)
		os.Exit(0)
	}
	if len(args) == 2 && args[1] == "--version" {
		fmt.Printf("fdvib %s\n", version)
		os.Exit(0)
	}
	if len(args) < 3 {
		usage(os.Stderr)
		os.Exit(2)
	}
	cmd := args[1]
	var err error
	switch {
	case cmd == "init" && len(args) == 3:
		cwd, getwdErr := os.Getwd()
		if getwdErr != nil {
			err = getwdErr
			break
		}
		err = initcmd.InitializeInput(args[2], cwd)
	case cmd == "-inp" && len(args) == 3:
		var s *settings.Settings
		s, err = settings.From(args[2], "")
		if err == nil {
			err = run.Calculate(s)
		}
	case cmd == "modes" && len(args) == 3:
		err = export.Modes(absolute(args[2]))
	case cmd == "shm" && len(args) == 3:
		err = export.Shm(absolute(args[2]))
	case cmd == "thermo" && len(args) == 5 && args[3] == "-inp":
		err = export.Thermo(absolute(args[2]), absolute(args[4]))
	default:
		usage(os.Stderr)
		os.Exit(2)
	}
	if err != nil {
		fmt.Fprintf(os.Stderr, "fdvib: error: %s\n", err)
		os.Exit(1)
	}
}

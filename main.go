package main

import (
	"fmt"
	"log"
	"os"
	"os/exec"
	"strings"

	"log/slog"

	"github.com/google/uuid"
)

func main() {
	programLevel := new(slog.LevelVar)
	logger := slog.New(slog.NewJSONHandler(os.Stderr, &slog.HandlerOptions{Level: programLevel}))
	slog.SetDefault(logger)

	integrationID := os.Getenv("INTEGRATION_ID")
	baseDir := os.Getenv("BASE_DIR")
	if integrationID == "" {
		id := uuid.New()
		integrationID = id.String()
	}
	if baseDir == "" {
		baseDir = "/mnt/efs"
	}

	logger.Info(integrationID)
	// create subdirectories
	err := os.Chdir(baseDir)
	if err != nil {
		logger.Error(err.Error())
		os.Exit(1)
	}

	// inputDir
	inputDir := fmt.Sprintf("%s/input/%s", baseDir, integrationID)
	err = os.MkdirAll(inputDir, 0755)
	if err != nil {
		logger.Error(err.Error())
		os.Exit(1)
	}

	// outputDir
	err = os.MkdirAll("output", 0777)
	if err != nil {
		logger.Error(err.Error())
		os.Exit(1)
	}
	err = os.Chown("output", 1000, 1000)
	if err != nil {
		logger.Error(err.Error())
		os.Exit(1)
	}

	outputDir := fmt.Sprintf("%s/output/%s", baseDir, integrationID)
	err = os.MkdirAll(outputDir, 0777)
	if err != nil {
		logger.Error(err.Error())
		os.Exit(1)
	}
	err = os.Chown(outputDir, 1000, 1000)
	if err != nil {
		logger.Error(err.Error())
		os.Exit(1)
	}

	log.Println("Starting pipeline")
	// run pipeline
	cmd := exec.Command("python3.9", "/service/main.py", inputDir, outputDir)
	cmd.Dir = "/service"
	var out strings.Builder
	var stderr strings.Builder
	cmd.Stdout = &out
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		logger.Error(err.Error(),
			slog.String("error", stderr.String()))
	}
	log.Println(out.String())

	logger.Info("Processing complete")
}

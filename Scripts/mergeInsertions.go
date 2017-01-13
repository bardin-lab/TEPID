package main

import (
	"encoding/csv"
	"flag"
	"fmt"
	"log"
	"os"
    "path/filepath"
	"strconv"
	"strings"
)

type Insertion struct {
	ins_chrom string
	ins_start int
	ins_end   int
	ref_chrom string
	ref_start int
	ref_end   int
	agi       []string
	accession []string
}

func readInsertion(path string) []Insertion {
	// Open BED file as CSV
	bedFile, err := os.Open(path)
	if err != nil {
		fmt.Println(err)
	}

	defer bedFile.Close()

	reader := csv.NewReader(bedFile)
	reader.Comma = '\t' // Use tab-delimited instead of comma <---- here!
	reader.FieldsPerRecord = -1

	bedData, err := reader.ReadAll()
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}

	var oneInsertion Insertion
	var allInsertions []Insertion

    // set sample identifier to input file - extension
    // i.e path = "test.bed" results in accession "test"
    basename := filepath.Base(path)
    ext := filepath.Ext(path)
    accession := strings.TrimRight(basename, ext)

	for _, each := range bedData {
		oneInsertion.ins_chrom = each[0]
		oneInsertion.ins_start, _ = strconv.Atoi(each[1]) // need to cast integer to string
		oneInsertion.ins_end, _ = strconv.Atoi(each[2])
		oneInsertion.ref_chrom = each[3]
		oneInsertion.ref_start, _ = strconv.Atoi(each[4])
		oneInsertion.ref_end, _ = strconv.Atoi(each[5])
		oneInsertion.agi = strings.Split(each[6], ",")
		oneInsertion.accession = strings.Fields(accession)
		allInsertions = append(allInsertions, oneInsertion)
	}
	return allInsertions
}

func sliceContains(sl []string, c string) bool {
	for _, s := range sl {
		if s == c {
			return true
		}
	}
	return false
}

func transposonsMatch(masterAgi, currentAgi []string) bool {
	for _, currentTE := range currentAgi {
		if sliceContains(masterAgi, currentTE) {
			return true
		}
	}
	return false
}

func overlap(masterInsertion, currentInsertion Insertion, extendLeftFlank, extendRightFlank int) bool {
	if masterInsertion.ins_start <= currentInsertion.ins_end+extendRightFlank {
        if currentInsertion.ins_start <= masterInsertion.ins_end+extendLeftFlank {
			return true
		}
	}
	return false
}

func canMergeExisting(master []Insertion, currentInsertion Insertion, extendLeftFlank, extendRightFlank int) bool {
	for i, masterInsertion := range master {
		if overlap(masterInsertion, currentInsertion, extendLeftFlank, extendRightFlank) {
			if transposonsMatch(masterInsertion.agi, currentInsertion.agi) {
				if masterInsertion.ins_chrom == currentInsertion.ins_chrom {
					//  adjust reference coords, taking coords of longest list of TEs. Why?
					if len(currentInsertion.agi) < len(masterInsertion.agi) {
						currentInsertion.ref_start, currentInsertion.ref_end = masterInsertion.ref_start, masterInsertion.ref_end
					}
					if masterInsertion.ref_start == currentInsertion.ref_start {
						if masterInsertion.ref_chrom == currentInsertion.ref_chrom {
							master[i] = mergeExisting(masterInsertion, currentInsertion)
							return true
						}
					}
				}
			}
		}
	}
	return false
}

func mergeExisting(masterInsertion, currentInsertion Insertion) Insertion {
	// merge accessions (=sample identifiers)
	for _, eachCurrentAccession := range currentInsertion.accession {
		if !sliceContains(masterInsertion.accession, eachCurrentAccession) {
			masterInsertion.accession = append(masterInsertion.accession, eachCurrentAccession)
		}
	}
	// merge agis (=inserted TE)
	for _, eachCurrentAgi := range currentInsertion.agi {
		if !sliceContains(masterInsertion.agi, eachCurrentAgi) {
			masterInsertion.agi = append(masterInsertion.agi, eachCurrentAgi)
		}
	}
	return masterInsertion
}

func mergeInsertions(master []Insertion, newInsertions []Insertion, extendLeftFlank, extendRightFlank int) []Insertion {
    // compares each new insertion with the complete master list.
    // could be refactored to sort and merge individiual files,
    // followed by iterating in pairs (or windows)
    // and evalutating whether the pairs (or windows)
    // can be joined. The current procedure does not guarantee
    // for the best possible merge and does not merge insertions
    // in the first sample.
	mergeable := 0
	for _, each := range newInsertions {
		canmerge := canMergeExisting(master, each, extendLeftFlank, extendRightFlank)
		if canmerge {
			mergeable = mergeable + 1
		} else {
			master = append(master, each)
		}
	}
	fmt.Printf("Found %d mergeable insertions\n", mergeable)
	return master
}

func writeResult(master []Insertion, result string) {
	file, err := os.Create(result)
	checkError("Cannot create file", err)
	defer file.Close()

	writer := csv.NewWriter(file)
	writer.Comma = '\t'

	line_template := "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s"

	for _, entry := range master {
		line := fmt.Sprintf(line_template, entry.ins_chrom, entry.ins_start, entry.ins_end, entry.ref_chrom, entry.ref_start, entry.ref_end, strings.Join(entry.agi, ","), strings.Join(entry.accession, ","))
		entry := strings.Split(line, "\t")
		err := writer.Write(entry)
		checkError("Cannot write to file", err)
	}

	defer writer.Flush()
}

func worker(samples []string, result string, extendLeftFlank, extendRightFlank int) {
	var Insertions [][]Insertion
	for _, each := range samples {
		insert := readInsertion(each)
		Insertions = append(Insertions, insert)
	}
	for _, each := range Insertions[1:] {
		Insertions[0] = mergeInsertions(Insertions[0], each, extendLeftFlank, extendRightFlank)
	}
	writeResult(Insertions[0], result)
}

func checkError(message string, err error) {
	if err != nil {
		log.Fatal(message, err)
	}
}

func (i *arrayFlags) String() string {
	return "my string representation"
}

func (i *arrayFlags) Set(value string) error {
	*i = append(*i, value)
	return nil
}

var inputFiles arrayFlags
var outFile string
var extendLeftFlank int
var extendRightFlank int

type arrayFlags []string

func main() {
	flag.Var(&inputFiles, "input", "Input files")
	flag.StringVar(&outFile, "outFile", "", "Output file")
	flag.IntVar(&extendLeftFlank, "extendLeftFlank", 0, "Extend flank by int nt")
	flag.IntVar(&extendRightFlank, "extendRightFlank", 100, "Extend flank by int nt")
	flag.Parse()
	if outFile == "" {
		fmt.Println("no outFile specified, exiting")
		return
	}
	if len(inputFiles) < 2 {
		fmt.Println("Need at least 2 input files to merge, exiting")
		return
	}
	worker(inputFiles, outFile, extendLeftFlank, extendRightFlank)
}

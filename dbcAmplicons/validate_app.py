# validate_app.py
#
import sys
import traceback
from dbcAmplicons import barcodeTable
from dbcAmplicons import primerTable
from dbcAmplicons import sampleTable


class validateApp:
    """
    Validate the sample, primer and barcode files
    """

    def __init__(self):
        self.verbose = False

    def validatePrimer(self, primerObject, debug=False):
        failed = False
        try:
            # Test completeness of primer file
            ptmp = primerObject.primers
            for p5i in primerObject.P5pair:
                try:
                    ptmp = list(set(ptmp) - set(primerObject.P5pair[p5i]))
                except ValueError:
                    pass
            if ptmp:
                sys.stderr.write("ERROR:[validate] %s pair ID(s) are missing from P5\n" % ','.join(ptmp))
                failed = True

            ptmp = primerObject.primers
            for p7i in primerObject.P7pair:
                try:
                    ptmp = list(set(ptmp) - set(primerObject.P7pair[p7i]))
                except ValueError:
                    pass
            if ptmp:
                sys.stderr.write("ERROR:[validate] %s pair ID(s) are missing from P7\n" % ','.join(ptmp))
                failed = True

            if failed:
                return 1
            return 0

        except (KeyboardInterrupt, SystemExit):
            self.clean()
            sys.stderr.write("%s unexpectedly terminated\n" % (__name__))
            return 1

        except:
            self.clean()
            if not debug:
                sys.stderr.write("A fatal error was encountered. trying turning on debug\n")
            if debug:
                sys.stderr.write("".join(traceback.format_exception(*sys.exc_info())))
            return 1

    def validateSample(self, barcodesObject, primerObject, samplesObject, debug=False):
        """
        Start preprocessing double barcoded Illumina sequencing run, perform
        """
        failed = False
        try:
            # Test validity of sample sheet
            for sampl_barcode in samplesObject.sampleTable:
                if sampl_barcode not in barcodesObject.IDS:
                    failed = True
                    sys.stderr.write("ERROR:[validate] barcode %s not found in barcode table\n" % sampl_barcode)
                for sampl_primer in samplesObject.sampleTable[sampl_barcode]:
                    if sampl_primer in ['-']:
                        continue
                    elif primerObject is None:
                        failed = True
                        sys.stderr.write("ERROR:[validate] primers found in sample sheet but no primer file provided\n")
                        raise
                    elif sampl_primer in ['*']:
                        continue
                    elif sampl_primer not in primerObject.primers:
                        failed = True
                        sys.stderr.write("ERROR:[validate] primer pair %s not found associated with barcode %s, sample %s in project %s\n"
                                % (sampl_primer, sampl_barcode,
                                samplesObject.sampleTable[sampl_barcode][sampl_primer][0],
                                samplesObject.sampleTable[sampl_barcode][sampl_primer][1]))
            if failed:
                return 1
            return 0

        except (KeyboardInterrupt, SystemExit):
            self.clean()
            sys.stderr.write("%s unexpectedly terminated\n" % (__name__))
            return 1

        except:
            self.clean()
            if not debug:
                sys.stderr.write("A fatal error was encountered. trying turning on debug\n")
            if debug:
                sys.stderr.write("".join(traceback.format_exception(*sys.exc_info())))
            return 1

    def start(self, barcodesFile, primerFile, samplesFile, verbose=True, debug=False):
        """
        Start preprocessing double barcoded Illumina sequencing run, perform
        """
        self.verbose = verbose
        try:
            # read in primer sequences
            bcTable = barcodeTable(barcodesFile)
            if self.verbose:
                sys.stdout.write("barcode table length: %s\n" % bcTable.getLength())

            # read in primer sequences if present
            prTable = primerTable(primerFile)
            if verbose:
                sys.stdout.write("primer table length P5 Primer (expanded) Sequences:%s, P7 Primer (expanded) Sequences:%s\n" % (len(prTable.getP5sequences()), len(prTable.getP7sequences())))
            res1 = self.validatePrimer(prTable, debug)

            # read in sample sheet
            sTable = sampleTable(samplesFile)
            if verbose:
                sys.stdout.write("sample table length: %s, and %s projects.\n" % (sTable.getSampleNumber(), len(sTable.getProjectList())))

            res2 = self.validateSample(bcTable, prTable, sTable, debug)
            if res1 == 0 and res2 == 0:
                sys.stderr.write("Validation confirmed, files are ok\n")
                return 0
            else:
                sys.stderr.write("Failed validation\n")
                return 1

        except (KeyboardInterrupt, SystemExit):
            self.clean()
            sys.stderr.write("%s unexpectedly terminated\n" % (__name__))
            return 1

        except:
            self.clean()
            if not debug:
                sys.stderr.write("A fatal error was encountered. trying turning on debug\n")
            if debug:
                sys.stderr.write("".join(traceback.format_exception(*sys.exc_info())))
            return 1

    def clean(self):
        if self.verbose:
            sys.stderr.write("Cleaning up.\n")

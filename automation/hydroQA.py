import os
import sys
import boto3
import argparse
import re
import logging
import time

"""

    Hydro QA Tool
    
    Non-standard Requirements: argparse, boto3

    usage: hydroQA.py [-h] [--verbose] [--dryrun] bucket imgdir
    
    positional arguments:
      bucket      Name of the program we are looking for.
      imgdir      Directory to load the 4 images.
    
    optional arguments:
      -h, --help  show this help message and exit
      --verbose   Give us more logs
      --dryrun    Don't touch S3


    Command line example: 
        c:\place\I\Pyt\Hydraulic-Modeling\automation\>python hydroQA.py sfr-champdata c:\imgs
        C:\Matt-SFR Files\Hydraulic Modeling\R Code to Build Input Files\R-Code\automation\>python hydroQA.py sfr-champdata c:\imgs
        C:\Hydro-QA\hydroQA.py --dryrun sfr-champdata c:\imgs


"""


def runHydroQA(bucket, imgdir, dryrun=False):
    # This is the main controller function that walks through the S3 bucket looking for QA results to
    # Work with and process.
    s3 = boto3.client('s3')

    # Loop over every object in the S3 bucket looking for project.rs.xml files in a Hydro/QA folder
    paginator = s3.get_paginator('list_objects')

    logging.info('Starting QA Process:::')

    for page in paginator.paginate(Bucket=bucket):
        if 'Contents' in page:
            for line in page['Contents']:
                # S3 doesn't give us a way to filter before the query so we just loop over everything and only
                # run the function when we've found a project inside a Hydro/QA folder
                if re.match(".*\/Hydro\/QA\/.*\/project.rs.xml$", line['Key']):
                    HydroQARun(bucket, imgdir, os.path.dirname(line['Key']), dryrun)

class HydroQARun():

    # Change these if you want different images downloaded:
    IMAGES = [
        'Depth.jpg',
        'Velocity.Magnitude.jpg',
        'Depth.Error.jpg',
        'Boundary_Conditions.jpg'
    ]
    def __init__(self, bucket, imgdir, basekey, dryrun):
        self.bucket = bucket
        self.imgdir = imgdir
        self.basekey = basekey
        self.dryrun = dryrun
        self.dryrunstr = "(dryrun)" if self.dryrun is True else ""

        logging.info("============================================================\n")
        logging.info("Result found: {}".format(self.basekey))

        # A little safety check to make sure we're not doing somethign stupid that will cause use to delete the entire
        # repo or something
        if not re.match(".*\/Hydro\/QA\/.*", self.basekey):
            logging.error("  Basekey was not correct: {}".format(self.basekey))
            return

        # Get our 4 images
        for imgname in HydroQARun.IMAGES:
            self.getLocalImage(imgname)

        # Do our ask.
        time.sleep(1)
        result = query_yes_no("  Images are downloaded. Is this Result Good?".format(basekey))

        # If we respond Y then move things from QA to Resultsw
        if result is True:
            logging.info("  QA Approved")
            self.cleanupDestFolder('Results')
            self.moveQAToFolder('Results')
        else:
            # If we respond "N" then move QA results to QA_Reject so they won't be picked up next time we run this script
            logging.info("  QA Rejected")
            self.cleanupDestFolder('QA_Reject')
            self.moveQAToFolder('QA_Reject')

        # Now clean things up for the next round of images
        for imgname in HydroQARun.IMAGES:
            self.deleteLocalImage(imgname)

    def cleanupDestFolder(self, foldername):
        """
        Clean up the destination folder in S3 so we don't have any old files polluting our new results
        :param foldername:
        :return:
        """

        s3 = boto3.client('s3')
        paginator = s3.get_paginator('list_objects')
        dst = self.basekey.replace('/Hydro/QA/', '/Hydro/{}/'.format(foldername))

        for page in paginator.paginate(Bucket=self.bucket, Prefix=dst):
            if 'Contents' in page:
                logging.info("Cleaning up existing S3 folder: {}".format(dst))
                for line in page['Contents']:
                    if self.dryrun is not True:
                        s3.delete_object(Bucket=self.bucket, Key=line['Key'])
                    logging.debug("    Deleted S3 Object: {} {}".format(line['Key'], self.dryrunstr))
            else:
                logging.info("    Nothing to clean up at: {}".format(dst))

    def moveQAToFolder(self, folder):
        """
        Move from the QA folder to whatever result folder we choose
        :param folder:
        :return:
        """

        s3 = boto3.client('s3')
        # Loop over every object in the S3 bucket looking for project.rs.xml files in a Hydro/QA folder
        paginator = s3.get_paginator('list_objects')
        for page in paginator.paginate(Bucket=self.bucket, Prefix=self.basekey):
            for line in page['Contents']:
                src = line['Key']
                dst = src.replace('/Hydro/QA', '/Hydro/{}'.format(folder))

                # Here's how move works in Boto
                if self.dryrun is not True:
                    s3.copy({"Bucket": self.bucket, "Key": src} , self.bucket, dst)
                    s3.delete_object(Bucket=self.bucket, Key=src)
                logging.debug("    Moved: {} ==> {} {}".format(src, dst, self.dryrunstr))


    def getLocalImage(self, imgname):
        """
        Download an image and put it in the the imgdir
        :param imgname:
        :return:
        """
        s3 = boto3.client('s3')
        localfile = os.path.join(self.imgdir, imgname)
        s3path = "/".join([self.basekey, imgname])
        # Delete the image if it already exists.
        if os.path.isfile(localfile):
            os.unlink(localfile)
        s3r = boto3.resource('s3')
        s3r.Bucket(self.bucket).download_file(s3path, localfile)
        # with open(localfile, 'w+') as f:
        #     f.write(s3.get_object(Bucket=self.bucket, Key=s3path)['Body'].read())
        logging.debug("  Downloaded Image: {}".format(imgname))

    def deleteLocalImage(self, imgname):
        """
        Delete the local images so we're sure what we're seeing is new.
        :return:
        """
        # Clean up images
        localfile = os.path.join(self.imgdir, imgname)
        if os.path.isfile(localfile):
            os.unlink(localfile)

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        print(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            print("Please respond with 'yes' or 'no' (or 'y' or 'n').\n")

def main():
    # This is just a really simple parser for command line options.

    parser = argparse.ArgumentParser()
    parser.add_argument('bucket',
                        help='Name of the program we are looking for.',
                        type=str)
    parser.add_argument('imgdir',
                        help='Directory to load the 4 images.',
                        type=str)
    parser.add_argument('--verbose',
                        action='store_true',
                        default=False,
                        help='Give us more logs')
    parser.add_argument('--dryrun',
                        action='store_true',
                        default=False,
                        help='Don\'t touch S3')
    args = parser.parse_args()

    try:

        # Make sure the output folder exists
        if not os.path.isdir(args.imgdir):
            os.makedirs(args.imgdir)

        # Set up logging
        logfile = os.path.join(args.imgdir, "hydroqa.log")

        if args.verbose:
            loglevel = logging.DEBUG
        else:
            loglevel = logging.INFO

        logging.getLogger('boto3').setLevel(logging.CRITICAL)
        logging.getLogger('botocore').setLevel(logging.CRITICAL)

        logFormatter = logging.Formatter("[%(levelname)-5.5s]  %(message)s")
        rootLogger = logging.getLogger()
        rootLogger.setLevel(loglevel)

        fileHandler = logging.FileHandler(logfile, mode='w')
        fileHandler.setFormatter(logFormatter)
        rootLogger.addHandler(fileHandler)

        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        rootLogger.addHandler(consoleHandler)

        runHydroQA(args.bucket, args.imgdir, args.dryrun)
    except AssertionError as e:
        print "Assertion Error", e
        sys.exit(0)
    except Exception as e:
        print 'Unexpected error: {0}'.format(sys.exc_info()[0]), e
        raise
        sys.exit(0)

"""
This handles the argument parsing and calls our main function
If we're not calling this from the command line then
"""
if __name__ == '__main__':
    main()

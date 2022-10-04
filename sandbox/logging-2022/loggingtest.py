import logging

format = "%(asctime)-15s [%(name)-12s] [%(levelname)-8s] %(message)s"
logging.basicConfig(format=format)
logging.addLevelName(25, 'PROGRESS')

logger = logging.getLogger('dtcc-builder')
logger.setLevel(logging.INFO)

logger.info('Computing building heights')
logger.warning('Height of building 15BE540D too small, skipping')
logger.log(25, '30.00%')
logger.info('Cleaning city model')
logger.info('Trimming 3D mesh')
logger.log(25, '30.00%')
logger.info('Extracting boundary of 3D mesh')
logger.info('Writing mesh to file CityMesh.vtu')
logger.log(25, '99.99%')

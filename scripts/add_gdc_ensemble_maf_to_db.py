import argparse
import csv
import gzip
import logging
from pathlib import Path
from typing import NamedTuple
from sqlalchemy import (
    create_engine, event,
    MetaData, Table, Column, Integer, Text,
)
from sqlalchemy.engine import Engine
from maf_utils import MAF

logger = logging.getLogger(__name__)
BATCH_SIZE = 5000


class GDCMutationCall(NamedTuple):
    sample: str
    chromosome: str
    start_position: int
    end_position: int
    strand: str
    reference_allele: str
    tumor_seq_allele1: str
    tumor_seq_allele2: str
    gdc_filter: str
    callers: str
    t_depth: int
    t_ref_count: int
    t_alt_count: int
    n_depth: int
    n_ref_count: int
    n_alt_count: int
    variant_type: str
    variant_classification: str
    hugo_symbol: str
    symbol: str
    hgnc_id: str
    gene: str
    transcript_id: str
    tsl: str
    biotype: str
    variant_class: str
    consequence: str
    one_consequence: str
    canonical: str
    hgvsc: str
    hgvsp: str
    hgvsp_short: str
    cdna_position: str
    cds_position: str
    protein_position: str
    amino_acids: str
    codons: str
    all_effects: str
    existing_variation: str
    cosmic: str

    @classmethod
    def create_cols(cls):
        """Create SQLAlchemy columns based on a named tuple class."""
        columns = [
            Column('id', Integer, primary_key=True, nullable=True)
        ]
        for field, f_type in cls._field_types.items():
            if f_type is int:
                col = Column(field, Integer)
            else:
                col = Column(field, Text)
            columns.append(col)
        return columns


    @classmethod
    def define_db_schema(cls, metadata, db_table_name):
        """Define the SQLite database schema."""
        Table(db_table_name, metadata, *cls.create_cols())


def read_maf(maf_pth, sample):
    """Read a MAF file."""
    def to_int(x):
        if x == '':
            return 0
        elif x == '.':
            return None
        else:
            return int(x)

    maf = MAF(maf_pth)
    for m in maf:
        yield GDCMutationCall(
            sample=sample, chromosome=m.chromosome,
            start_position=int(m.start_position), end_position=int(m.end_position),
            strand=m.strand,
            reference_allele=m.reference_allele,
            tumor_seq_allele1=m.tumor_seq_allele1,
            tumor_seq_allele2=m.tumor_seq_allele2,
            gdc_filter=m.gdc_filter, callers=m.callers,
            t_depth=to_int(m.t_depth), t_ref_count=to_int(m.t_ref_count),
            t_alt_count=to_int(m.t_alt_count),
            n_depth=to_int(m.n_depth), n_ref_count=to_int(m.n_ref_count),
            n_alt_count=to_int(m.n_alt_count),
            variant_type=m.variant_type, variant_classification=m.variant_classification,
            hugo_symbol=m.hugo_symbol, symbol=m.symbol,
            gene=m.gene, hgnc_id=m.hgnc_id,
            transcript_id=m.transcript_id, tsl=m.tsl,
            biotype=m.biotype, variant_class=m.variant_class,
            consequence=m.consequence, one_consequence=m.one_consequence,
            canonical=m.canonical,
            hgvsc=m.hgvsc, hgvsp=m.hgvsp, hgvsp_short=m.hgvsp_short,
            cdna_position=m.cdna_position, cds_position=m.cds_position, protein_position=m.protein_position,
            amino_acids=m.amino_acids, codons=m.codons,
            all_effects=m.all_effects,
            cosmic=m.cosmic, existing_variation=m.existing_variation
        )


def read_all_mafs(manifest_pth):
    """Sequentially read all MAFs."""
    with gzip.open(manifest_pth, 'rt') as f:
        manifest_reader = csv.DictReader(f, dialect='excel-tab')
        mafs = [
            (str(Path('external_data/cptac3_aliquot_ensemble_masked_mafs', row['maf_filename'])), row['SUBJECT_ID'])
            for row in manifest_reader
        ]

    logger.info(f'Load total {len(mafs)} MAFs')

    num_muts = 0
    for maf_pth, sample in mafs:
        logger.info(f'... loading {maf_pth}')
        for i, mut in enumerate(read_maf(maf_pth, sample), 1):
            if i % 50000 == 0:
                logger.info(f'... processed {i:,d} mutation calls')
            yield mut
        logger.info(f'... loaded {i:,d} mutation calls')
        num_muts += i
    logger.info(f'Loaded total {num_muts:,d} mutation calls')


def load_muts_to_db(conn, metadata, muts, db_table_name):
    """Load mutations to a given table.

    muts is an iterable of MAFMutationCall objects.
    """
    ins = metadata.tables[db_table_name].insert()
    # Insert mutation calls by batches
    ins_batch = []
    for mut in muts:
        # Add new record into batch
        ins_batch.append(mut._asdict())

        # Commit the batch when a batch is full
        if len(ins_batch) >= BATCH_SIZE:
            with conn.begin():
                conn.execute(ins, ins_batch)
            ins_batch = []

    # Load the last batch
    if ins_batch:
        with conn.begin():
            conn.execute(ins, ins_batch)


@event.listens_for(Engine, "connect")
def set_sqlite_pragma(dbapi_connection, connection_record):
    cursor = dbapi_connection.cursor()
    cursor.execute("PRAGMA cache_size=-4000000")
    cursor.execute("PRAGMA temp_store=MEMORY")
    cursor.execute("PRAGMA journal_mode=MEMORY")
    cursor.close()


def setup_cli():
    # Setup console logging
    console = logging.StreamHandler()
    all_loggers = logging.getLogger()
    all_loggers.setLevel(logging.INFO)
    all_loggers.addHandler(console)
    log_fmt = '[%(asctime)s][%(levelname)-7s] %(message)s'
    log_formatter = logging.Formatter(log_fmt, '%Y-%m-%d %H:%M:%S')
    console.setFormatter(log_formatter)

    parser = argparse.ArgumentParser(
        description='Load GDC ensemble MAF into a SQLite database'
    )
    parser.add_argument('db_pth', help="Path to the SQLite database")
    parser.add_argument('table', help="Table name")
    parser.add_argument('manifest_pth', help="Path to the manifest file")
    return parser


def main(db_pth, manifest_pth, db_table_name):
    logger.info(f'Add new table {db_table_name} to the SQLite database at {db_pth}')
    metadata = MetaData()
    db_engine = create_engine(f'sqlite:///{db_pth}')
    GDCMutationCall.define_db_schema(metadata, db_table_name)
    metadata.create_all(db_engine)
    conn = db_engine.connect()

    all_muts = read_all_mafs(manifest_pth)
    load_muts_to_db(conn, metadata, all_muts, db_table_name)
    logger.info('Complete')


if __name__ == '__main__':
    parser = setup_cli()
    args = parser.parse_args()
    main(args.db_pth, args.manifest_pth, args.table)


    # --------------------------------------------------------
    # PROCESS EACH REPLICON
    # --------------------------------------------------------

    for accession in common_accessions:

        print("\n" + "=" * 65)
        print("PROCESSING:", accession)
        print("=" * 65)

        original_record = original_records[accession]
        rotated_record = rotated_records[accession]

        original_temp = None
        rotated_temp = None

        try:

            # ------------------------------------------------
            # CREATE TEMPORARY SINGLE-RECORD FASTA FILES
            # ------------------------------------------------

            with tempfile.NamedTemporaryFile(
                mode="w",
                suffix=".fna",
                prefix=f".{accession}_original_",
                dir=output_dir,
                delete=False
            ) as temp:

                original_temp = temp.name

            with tempfile.NamedTemporaryFile(
                mode="w",
                suffix=".fna",
                prefix=f".{accession}_rotated_",
                dir=output_dir,
                delete=False
            ) as temp:

                rotated_temp = temp.name

            # Write original record
            SeqIO.write(
                original_record,
                original_temp,
                "fasta"
            )

            # Write corresponding rotated record
            SeqIO.write(
                rotated_record,
                rotated_temp,
                "fasta"
            )

            # ------------------------------------------------
            # GENERATE DOT PLOT
            # ------------------------------------------------

            output_png = plot_dotplot(
                os.path.basename(original_temp),
                os.path.basename(rotated_temp),
                report_file,
                output_dir,
                plots_dir
            )

            print("PNG:", output_png)

            print(
                "SVG:",
                output_png.replace(".png", ".svg")
            )

            successful += 1

        except Exception as e:

            print(
                f"ERROR processing {accession}: {e}"
            )

            failed += 1

        finally:

            # --------------------------------------------
            # REMOVE TEMPORARY FASTA FILES
            # --------------------------------------------

            if original_temp and os.path.exists(original_temp):
                os.remove(original_temp)

            if rotated_temp and os.path.exists(rotated_temp):
                os.remove(rotated_temp)

    # --------------------------------------------------------
    # FINAL SUMMARY
    # --------------------------------------------------------

    print("\n" + "=" * 65)
    print("DOT PLOT GENERATION COMPLETED")
    print("=" * 65)

    print("Total original replicons:", len(original_records))
    print("Matching replicons:", len(common_accessions))
    print("Successfully generated:", successful)
    print("Failed:", failed)

    print("\nOutput directory:", plots_dir)

